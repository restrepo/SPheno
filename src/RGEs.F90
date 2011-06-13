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
! new variables
!------------------------------------------------
 Integer, Private, Save :: loop_ord
 Logical, Private, Save :: Nu_R, Phi, S_f, LH, LLE, UDD, l_mu, l_lin, l_qua
!------------------------------------------------
! needed for programs by Florian Staub
!------------------------------------------------
# ifdef SARAH
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

# endif SARAH
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

  diagonal(1,1) = 3._dp * TraceY(2) + TraceY(1)     &
              & + c1_1(1,1) * gauge2(1) + c1_1(1,2) * gauge2(2)
  sume1 = 3._dp * aYeYe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do

  betaYe1 = Matmul2(Ye,sume1,OnlyDiagonal)

  diagonal(2,1) = 3._dp * TraceY(2) + TraceY(1)              &
    &  + c1_1(2,1) * gauge2(1) + c1_1(2,2) * gauge2(2) + c1_1(2,3) * gauge2(3)
  sumd1  = 3._dp * aYdYd + aYuYu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do

  betaYd1 = Matmul2(Yd,sumd1,OnlyDiagonal)

  diagonal(3,1) = 3._dp * TraceY(3)              &
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
!   aYdYdaYuYu = Matmul2(aYdYd,aYuYu,OnlyDiagonal)
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

  !------------------------------------
  ! beta functions for A-parameters
  !-----------------------------------
  Call Adjungate(Ad,aAd)
  Call Adjungate(Ae,aAe)
  Call Adjungate(Au,aAu)

  aAdAd = MatMul2(aAd,Ad,OnlyDiagonal)
  aAeAe = MatMul2(aAe,Ae,OnlyDiagonal)
  aAuAu = MatMul2(aAu,Au,OnlyDiagonal)

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

  aYdAd = MatMul2(aYd,Ad,OnlyDiagonal)
  aYeAe = MatMul2(aYe,Ae,OnlyDiagonal)
  aYuAu = MatMul2(aYu,Au,OnlyDiagonal)

  TraceaYA(1) = cTrace(aYeAe) 
  TraceaYA(2) = cTrace(aYdAd) 
  TraceaYA(3) = cTrace(aYuAu) 

  g2Mi = gauge2 * Mhlf
  !--------------
  ! A_e
  !--------------
  sume1 = sume1 + 2._dp * aYeYe
  betaAe1 = MatMul2(Ae,sume1,OnlyDiagonal)

  diagonal(1,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1)  &
                &         - c1_1(1,1) * g2Mi(1) - c1_1(1,2) * g2Mi(2)    )
  sume1 = 4._dp * aYeAe
  Do i1=1,3
   sume1(i1,i1) = sume1(i1,i1) + diagonal(1,1)
  End Do 
  betaAe1 = betaAe1 + MatMul2(Ye,sume1,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
  sumd1 = sumd1 + 2._dp * aYdYd
  betaAd1 = MatMul2(Ad,sumd1,OnlyDiagonal)
  
  diagonal(2,1) = 2._dp * ( 3._dp * TraceaYA(2) + TraceaYA(1) &
                &         - c1_1(2,1) * g2Mi(1) - c1_1(2,2) * g2Mi(2)   &
                &         - c1_1(2,3) * g2Mi(3) )
  sumd1 = 4._dp * aYdAd + 2._dp * aYuAu
  Do i1=1,3
   sumd1(i1,i1) = sumd1(i1,i1) + diagonal(2,1)
  End Do
  betaAd1 = betaAd1 + MatMul2(Yd,sumd1,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
  sumu1 = sumu1 + 2._dp * aYuYu
  betaAu1 = MatMul2(Au,sumu1,OnlyDiagonal)
  
  diagonal(3,1) = 2._dp * ( 3._dp * TraceaYA(3)              &
                &         - c1_1(3,1) * g2Mi(1) - c1_1(3,2) * g2Mi(2)   &
                &         - c1_1(3,3) * g2Mi(3) )
  sumu1 = 2._dp * aYdAd + 4._dp * aYuAu
  Do i1=1,3
   sumu1(i1,i1) = sumu1(i1,i1) + diagonal(3,1)
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
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)   &
       & - 6._dp * gauge2(2) + 1.2_dp * gauge2(1)
   sume2 = sume2 - 2._dp * aYeYeaYeYe - hd(1) * aYeYe
   betaAe2 = MatMul2(Ae,sume2,OnlyDiagonal)
    
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
   betaAe2 = betaAe2 + MatMul2(Ye,sume2,OnlyDiagonal)

  !--------------
  ! A_d
  !--------------
   hd(1) = 6._dp * TraceY(2) + 2._dp * TraceY(1)    &
       & - 6._dp * gauge2(2) - 0.4_dp * gauge2(1)
   sumd2 = sumd2 - 2._dp * ( aYdYdaYdYd + aYuYuaYdYd ) - hd(1) * aYdYd
   betaAd2 = MatMul2(Ad,sumd2,OnlyDiagonal)
    
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
   betaAd2 = betaAd2 + MatMul2(Yd,sumd2,OnlyDiagonal)

  !--------------
  ! A_u
  !--------------
   hd(1) = 6._dp * ( TraceY(3) - gauge2(2) ) + 0.4_dp * gauge2(1)
   sumu2 = sumu2 - 2._dp * ( aYuYuaYuYu + aYdYdaYuYu ) - hd(1) * aYuYu
   betaAu2 = MatMul2(Au,sumu2,OnlyDiagonal)
    
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
   YeaYe = MatMul2(Ye,aYe,OnlyDiagonal)
   YuaYu = MatMul2(Yu,aYu,OnlyDiagonal)

   MdYdaYd = MatMul2(Md,YdaYd,OnlyDiagonal)
   MeYeaYe = MatMul2(Me,YeaYe,OnlyDiagonal)
   MlaYeYe = MatMul2(Ml,aYeYe,OnlyDiagonal)
   MqaYdYd = MatMul2(Mq,aYdYd,OnlyDiagonal)
   MqaYuYu = MatMul2(Mq,aYuYu,OnlyDiagonal)
   MuYuaYu = MatMul2(Mu,YuaYu,OnlyDiagonal)

   Call Adjungate(MdYdaYd, YdaYdMd) ! YdaYdMd = MatMul2(YdaYd,Md,OnlyDiagonal)
   Call Adjungate(MeYeaYe, YeaYeMe) ! YeaYeMe = MatMul2(YeaYe,Me,OnlyDiagonal)
   Call Adjungate(MlaYeYe, aYeYeMl) ! aYeYeMl = MatMul2(aYeYe,Ml,OnlyDiagonal)
   Call Adjungate(MqaYdYd, aYdYdMq) ! aYdYdMq = MatMul2(aYdYd,Mq,OnlyDiagonal)
   Call Adjungate(MqaYuYu, aYuYuMq) ! aYuYuMq = MatMul2(aYuYu,Mq,OnlyDiagonal)
   Call Adjungate(MuYuaYu, YuaYuMu) ! YuaYuMu = MatMul2(YuaYu,Mu,OnlyDiagonal)

   aYdMdYd = MatMul3(aYd,Md,Yd,OnlyDiagonal)
   aYeMeYe = MatMul3(aYe,Me,Ye,OnlyDiagonal)
   aYuMuYu = MatMul3(aYu,Mu,Yu,OnlyDiagonal)
   YdMqaYd = MatMul3(Yd,Mq,aYd,OnlyDiagonal)
   YeMlaYe = MatMul3(Ye,Ml,aYe,OnlyDiagonal)
   YuMqaYu = MatMul3(Yu,Mq,aYu,OnlyDiagonal)

   AdaAd = MatMul2(Ad,aAd,OnlyDiagonal)
   AeaAe = MatMul2(Ae,aAe,OnlyDiagonal)
   AuaAu = MatMul2(Au,aAu,OnlyDiagonal)

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

    AdaYd = MatMul2(Ad,aYd,OnlyDiagonal)
    AeaYe = MatMul2(Ae,aYe,OnlyDiagonal)
    AuaYu = MatMul2(Au,aYu,OnlyDiagonal)

    aAdYd = MatMul2(aAd,Yd,OnlyDiagonal)
    aAeYe = MatMul2(aAe,Ye,OnlyDiagonal)
    aAuYu = MatMul2(aAu,Yu,OnlyDiagonal)

    Call Adjungate(AdaYd,YdaAd) ! YdaAd = MatMul2(Yd,aAd,OnlyDiagonal)
    Call Adjungate(AeaYe,YeaAe) ! YeaAe = MatMul2(Ye,aAe,OnlyDiagonal)
    Call Adjungate(AuaYu,YuaAu) ! YuaAu = MatMul2(Yu,aAu,OnlyDiagonal)

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
     & - 3._dp * (MH(2)*TraceY(3) - MH(1) * TraceY(2) )        &
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

    TraceB(1) = cTrace( 3._dp * (Matmul2(AuaYu,YuaYu,OnlyDiagonal)     &
              &                 + Matmul2(AdaYd,YdaYd,OnlyDiagonal) )  &
              &       + MatMul2(AeaYe,YeaYe,OnlyDiagonal)              &
              &       + MatMul2(aYuAu,aYdYd,OnlyDiagonal)              &
              &       + MatMul2(aYdAd,aYuYu,OnlyDiagonal) ) 
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

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, YTaYT, betaYT1, sumT1     &
      & , AT, aAT, aATAT, ATaAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl &
      & , DYT, DAT
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
   YTaYT = MatMul2(YT,aYT,OnlyDiagonal)
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
   ATaAT = MatMul2(AT,aAT,OnlyDiagonal)
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

  Complex(dp), Dimension(3,3) :: YT, aYT, aYTYT, YTaYT, betaYT1, sumT1     &
      & , AT, aAT, aATAT, ATaAT, betaAT1, aYTAT, MlaYTYT, aYTMlYT, aYTYTMl &
      & , DYT, DAT
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
   YTaYT = MatMul2(YT,aYT,OnlyDiagonal)
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
   ATaAT = MatMul2(AT,aAT,OnlyDiagonal)
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
  Complex(dp), Dimension(3,3) :: A_u, A_d, A_N, AdaYd, aYdAd, AnaYn, aYnAn &
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
  AnaYn = Matmul2(A_n, aY_n,OnlyDiagonal)
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


#ifdef SARAH

 Subroutine GToParameters111(g,g1,g2,g3,Yu,Yd,Ye,Yb3,Yw3,Yx3)

 Implicit None 
  Real(dp), Intent(in) :: g(111) 
  Real(dp),Intent(out) :: g1,g2,g3

  Complex(dp), Intent(out), Dimension(3,3) :: Yu, Yd, Ye, Yb3, Yw3, Yx3

  Integer i1, i2, i3, i4, SumI 
 
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

  Integer i1, i2, i3, i4, SumI 
 
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

Subroutine rge111(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2,i3,i4 
Integer :: j1,j2,j3,j4,j5,j6,j7 
Real(dp) :: q 
Real(dp) :: g1,betag11,betag12,Dg1,g2,betag21,betag22,Dg2,g3,betag31,betag32,Dg3
Complex(dp) :: Yu(3,3),betaYu1(3,3),betaYu2(3,3),DYu(3,3),adjYu(3,3),Yd(3,3)   & 
& ,betaYd1(3,3),betaYd2(3,3),DYd(3,3),adjYd(3,3),Ye(3,3),betaYe1(3,3),betaYe2(3,3)       & 
& ,DYe(3,3),adjYe(3,3),Yb3(3,3),betaYb31(3,3),betaYb32(3,3),DYb3(3,3),adjYb3(3,3)        & 
& ,Yw3(3,3),betaYw31(3,3),betaYw32(3,3),DYw3(3,3),adjYw3(3,3),Yx3(3,3),betaYx31(3,3)     & 
& ,betaYx32(3,3),DYx3(3,3),adjYx3(3,3)
Complex(dp) :: Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3),Yx3adjYx3(3,3)   & 
& ,adjYb3Yb3(3,3),adjYdYd(3,3),adjYeYe(3,3),adjYuYu(3,3),adjYw3Yw3(3,3),adjYx3Yx3(3,3)   & 
& ,CYdTpYd(3,3),CYx3Yd(3,3),Yb3adjYb3Yb3(3,3),Yb3adjYeYe(3,3),Yb3adjYw3Yw3(3,3)          & 
& ,YdadjYdYd(3,3),YdadjYuYu(3,3),YeadjYb3Yb3(3,3),YeadjYeYe(3,3),YeadjYw3Yw3(3,3)        & 
& ,YuadjYdYd(3,3),YuadjYuYu(3,3),Yw3adjYb3Yb3(3,3),Yw3adjYeYe(3,3),Yw3adjYw3Yw3(3,3)     & 
& ,Yx3adjYx3Yx3(3,3),Yx3CYdTpYd(3,3),TpYx3CYx3Yd(3,3)

Complex(dp) :: Yb3adjYe(3,3),Yb3adjYw3(3,3),YdadjYu(3,3),YeadjYb3(3,3),YeadjYw3(3,3),YuadjYd(3,3)    & 
& ,Yw3adjYb3(3,3),Yw3adjYe(3,3),CYuTpYd(3,3),TpYx3CYx3(3,3),adjYb3Yb3adjYb3(3,3)         & 
& ,adjYb3Yb3adjYe(3,3),adjYb3Yb3adjYw3(3,3),adjYdYdadjYd(3,3),adjYdYdadjYu(3,3)          & 
& ,adjYdTpYx3CYx3(3,3),adjYeYeadjYb3(3,3),adjYeYeadjYe(3,3),adjYeYeadjYw3(3,3)           & 
& ,adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYw3Yw3adjYb3(3,3),adjYw3Yw3adjYe(3,3)          & 
& ,adjYw3Yw3adjYw3(3,3),adjYx3Yx3adjYx3(3,3),CYx3YdadjYd(3,3),TpYdadjYx3Yx3(3,3)         & 
& ,TpYdCYdTpYd(3,3),TpYuCYuTpYd(3,3),Yb3adjYb3Yb3adjYb3(3,3),Yb3adjYeYeadjYb3(3,3)       & 
& ,Yb3adjYw3Yw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYuYuadjYd(3,3) & 
& ,YeadjYb3Yb3adjYe(3,3),YeadjYeYeadjYe(3,3),YeadjYw3Yw3adjYe(3,3),YuadjYdYdadjYu(3,3)   & 
& ,YuadjYuYuadjYu(3,3),Yw3adjYb3Yb3adjYw3(3,3),Yw3adjYeYeadjYw3(3,3),Yw3adjYw3Yw3adjYw3(3,3)& 
& ,Yx3adjYx3Yx3adjYx3(3,3),adjYb3Yb3adjYb3Yb3(3,3),adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYw3Yw3(3,3)& 
& ,adjYdYdadjYdYd(3,3),adjYdYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYeYeadjYb3Yb3(3,3)   & 
& ,adjYeYeadjYeYe(3,3),adjYeYeadjYw3Yw3(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYuYu(3,3)     & 
& ,adjYw3Yw3adjYb3Yb3(3,3),adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYw3Yw3(3,3),adjYx3Yx3adjYx3Yx3(3,3)& 
& ,CYdTpYdadjYx3Yx3(3,3),CYdTpYdCYdTpYd(3,3),CYdTpYuCYuTpYd(3,3),CYx3TpYx3CYx3Yd(3,3)    & 
& ,TpYx3CYx3YdadjYd(3,3),Yb3adjYb3Yb3adjYb3Yb3(3,3),Yb3adjYb3Yb3adjYw3Yw3(3,3)           & 
& ,Yb3adjYeYeadjYb3Yb3(3,3),Yb3adjYeYeadjYeYe(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3)            & 
& ,Yb3adjYw3Yw3adjYw3Yw3(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3)              & 
& ,YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYuYu(3,3),YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYeYe(3,3)& 
& ,YeadjYb3Yb3adjYw3Yw3(3,3),YeadjYeYeadjYeYe(3,3),YeadjYw3Yw3adjYb3Yb3(3,3)             & 
& ,YeadjYw3Yw3adjYeYe(3,3),YeadjYw3Yw3adjYw3Yw3(3,3),YuadjYdYdadjYdYd(3,3)               & 
& ,YuadjYdYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),YuadjYuYuadjYuYu(3,3),Yw3adjYb3Yb3adjYb3Yb3(3,3)& 
& ,Yw3adjYb3Yb3adjYw3Yw3(3,3),Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYw3Yw3(3,3)            & 
& ,Yw3adjYw3Yw3adjYb3Yb3(3,3),Yw3adjYw3Yw3adjYw3Yw3(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3)      & 
& ,Yx3CYdTpYdadjYx3Yx3(3,3),Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYuCYuTpYd(3,3),               & 
& TpYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: TrYb3adjYb3,TrYdadjYd,TrYeadjYe,TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3

Complex(dp) :: TrYb3adjYb3Yb3adjYb3,TrYb3adjYeYeadjYb3,TrYb3adjYw3Yw3adjYb3,TrYdadjYdYdadjYd,        & 
& TrYdadjYdTpYx3CYx3,TrYdadjYuYuadjYd,TrYeadjYb3Yb3adjYe,TrYeadjYeYeadjYe,               & 
& TrYeadjYw3Yw3adjYe,TrYuadjYdYdadjYu,TrYuadjYuYuadjYu,TrYw3adjYb3Yb3adjYw3,             & 
& TrYw3adjYeYeadjYw3,TrYw3adjYw3Yw3adjYw3,TrYx3adjYx3Yx3adjYx3,TrTpYx3CYx3YdadjYd

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g1p5,g2p4,g2p5,g3p4,g3p5

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
 Yb3adjYe = Matmul2(Yb3,adjYe,OnlyDiagonal) 
 Yb3adjYw3 = Matmul2(Yb3,adjYw3,OnlyDiagonal) 
 YdadjYu = Matmul2(Yd,adjYu,OnlyDiagonal) 
 YeadjYb3 = Matmul2(Ye,adjYb3,OnlyDiagonal) 
 YeadjYw3 = Matmul2(Ye,adjYw3,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 Yw3adjYb3 = Matmul2(Yw3,adjYb3,OnlyDiagonal) 
 Yw3adjYe = Matmul2(Yw3,adjYe,OnlyDiagonal) 
 CYuTpYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 TpYx3CYx3 = Matmul2(Transpose(Yx3),Conjg(Yx3),OnlyDiagonal) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 adjYb3Yb3adjYb3 = Matmul2(adjYb3,Yb3adjYb3,OnlyDiagonal) 
 adjYb3Yb3adjYe = Matmul2(adjYb3,Yb3adjYe,OnlyDiagonal) 
 adjYb3Yb3adjYw3 = Matmul2(adjYb3,Yb3adjYw3,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdYdadjYu = Matmul2(adjYd,YdadjYu,OnlyDiagonal) 
 adjYdTpYx3CYx3 = Matmul2(adjYd,TpYx3CYx3,OnlyDiagonal) 
 adjYeYeadjYb3 = Matmul2(adjYe,YeadjYb3,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYeYeadjYw3 = Matmul2(adjYe,YeadjYw3,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYw3Yw3adjYb3 = Matmul2(adjYw3,Yw3adjYb3,OnlyDiagonal) 
 adjYw3Yw3adjYe = Matmul2(adjYw3,Yw3adjYe,OnlyDiagonal) 
 adjYw3Yw3adjYw3 = Matmul2(adjYw3,Yw3adjYw3,OnlyDiagonal) 
 adjYx3Yx3adjYx3 = Matmul2(adjYx3,Yx3adjYx3,OnlyDiagonal) 
 CYx3YdadjYd = Matmul2(Conjg(Yx3),YdadjYd,OnlyDiagonal) 
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
 YeadjYb3Yb3adjYe = Matmul2(Ye,adjYb3Yb3adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYb3Yb3adjYe(i2,i2) =  Real(YeadjYb3Yb3adjYe(i2,i2),dp) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYw3Yw3adjYe = Matmul2(Ye,adjYw3Yw3adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 Yw3adjYb3Yb3adjYw3 = Matmul2(Yw3,adjYb3Yb3adjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYb3Yb3adjYw3(i2,i2) =  Real(Yw3adjYb3Yb3adjYw3(i2,i2),dp) 
 Yw3adjYeYeadjYw3 = Matmul2(Yw3,adjYeYeadjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYeYeadjYw3(i2,i2) =  Real(Yw3adjYeYeadjYw3(i2,i2),dp) 
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
 TpYx3CYx3YdadjYd = Matmul2(Transpose(Yx3),CYx3YdadjYd,OnlyDiagonal) 
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
 TrYeadjYb3Yb3adjYe = cTrace(YeadjYb3Yb3adjYe) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYuadjYdYdadjYu = cTrace(YuadjYdYdadjYu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYw3adjYb3Yb3adjYw3 = cTrace(Yw3adjYb3Yb3adjYw3) 
 TrYw3adjYeYeadjYw3 = cTrace(Yw3adjYeYeadjYw3) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 TrTpYx3CYx3YdadjYd = cTrace(TpYx3CYx3YdadjYd) 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 33._dp*g1p3/5._dp + 5._dp*g1p3*NGHx3/2._dp + 5._dp*g1p3*NGHxb3/2._dp

 
 
If (TwoLoopRGE) Then 
betag12 = 199._dp*g1p5/25._dp + 27._dp*g1p3*g2p2/5._dp + 88._dp*g1p3*g3p2/5._dp -               & 
&  9._dp*g1p3*TrYb3adjYb3/25._dp - 14._dp*g1p3*TrYdadjYd/5._dp - 18._dp*g1p3*TrYeadjYe/5._dp -& 
&  26._dp*g1p3*TrYuadjYu/5._dp - 12._dp*g1p3*TrYw3adjYw3/5._dp - 38._dp*g1p3*TrYx3adjYx3/5._dp +& 
&  25._dp*g1p5*NGHx3/6._dp + 15._dp*g1p3*g2p2*NGHx3/2._dp +& 
&  40._dp*g1p3*g3p2*NGHx3/3._dp + 25._dp*g1p5*NGHxb3/6._dp +& 
&  15._dp*g1p3*g2p2*NGHxb3/2._dp + 40._dp*g1p3*g3p2*NGHxb3/3._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = g2p3 + 2._dp*g2p3*NGHw3 + 3._dp*g2p3*NGHx3/2._dp +& 
&  3._dp*g2p3*NGHxb3/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = 9._dp*g1p2*g2p3/5._dp + 25._dp*g2p5 + 24._dp*g2p3*g3p2 - 3._dp*g2p3*TrYb3adjYb3/5._dp -& 
&  6._dp*g2p3*TrYdadjYd - 2._dp*g2p3*TrYeadjYe - 6._dp*g2p3*TrYuadjYu - 28._dp*g2p3*TrYw3adjYw3/3._dp -& 
&  6._dp*g2p3*TrYx3adjYx3 + 24._dp*g2p5*NGHw3 + 5._dp*g1p2*g2p3*NGHx3/2._dp +& 
&  21._dp*g2p5*NGHx3/2._dp + 8._dp*g2p3*g3p2*NGHx3 +       & 
&  5._dp*g1p2*g2p3*NGHxb3/2._dp + 21._dp*g2p5*NGHxb3/2._dp +& 
&  8._dp*g2p3*g3p2*NGHxb3

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = -3._dp*g3p3 + 3._dp*g3p3*NGHg3 + g3p3*NGHx3 +& 
&  g3p3*NGHxb3

 
 
If (TwoLoopRGE) Then 
betag32 = 11._dp*g1p2*g3p3/5._dp + 9._dp*g2p2*g3p3 + 14._dp*g3p5 - 4._dp*g3p3*TrYdadjYd -       & 
&  4._dp*g3p3*TrYuadjYu - 4._dp*g3p3*TrYx3adjYx3 + 54._dp*g3p5*NGHg3 +    & 
&  5._dp*g1p2*g3p3*NGHx3/3._dp + 3._dp*g2p2*g3p3*NGHx3 +   & 
&  34._dp*g3p5*NGHx3/3._dp + 5._dp*g1p2*g3p3*NGHxb3/3._dp +& 
&  3._dp*g2p2*g3p3*NGHxb3 + 34._dp*g3p5*NGHxb3/3._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = -13._dp*g1p2*Yu/15._dp - 3._dp*g2p2*Yu - 16._dp*g3p2*Yu/3._dp +            & 
&  3._dp*TrYb3adjYb3*Yu/10._dp + 3._dp*TrYuadjYu*Yu + 3._dp*TrYw3adjYw3*Yu/2._dp +       & 
&  3._dp*TrYx3adjYx3*Yu + YuadjYdYd + 3._dp*YuadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYu2 = 2743._dp*g1p4*Yu/450._dp + g1p2*g2p2*Yu + 15._dp*g2p4*Yu/2._dp + 136._dp*g1p2*g3p2*Yu/45._dp +& 
&  8._dp*g2p2*g3p2*Yu - 16._dp*g3p4*Yu/9._dp - 27._dp*TrYb3adjYb3Yb3adjYb3*Yu/100._dp -  & 
&  3._dp*TrYb3adjYeYeadjYb3*Yu/10._dp - 3._dp*TrYb3adjYw3Yw3adjYb3*Yu/4._dp -            & 
&  6._dp*TrYdadjYdTpYx3CYx3*Yu - 3._dp*TrYuadjYdYdadjYu*Yu + 4._dp*g1p2*TrYuadjYu*Yu/5._dp +& 
&  16._dp*g3p2*TrYuadjYu*Yu - 9._dp*TrYuadjYuYuadjYu*Yu - 27._dp*TrYw3adjYb3Yb3adjYw3*Yu/40._dp -& 
&  3._dp*TrYw3adjYeYeadjYw3*Yu/2._dp + 6._dp*g2p2*TrYw3adjYw3*Yu - 15._dp*TrYw3adjYw3Yw3adjYw3*Yu/4._dp +& 
&  2._dp*g1p2*TrYx3adjYx3*Yu + 16._dp*g3p2*TrYx3adjYx3*Yu - 9._dp*TrYx3adjYx3Yx3adjYx3*Yu -& 
&  2._dp*YuadjYdTpYx3CYx3Yd + 2._dp*g1p2*YuadjYdYd/5._dp - 3._dp*TrYdadjYd*YuadjYdYd -   & 
&  TrYeadjYe*YuadjYdYd - 2._dp*YuadjYdYdadjYdYd - 2._dp*YuadjYdYdadjYuYu +               & 
&  2._dp*g1p2*YuadjYuYu/5._dp + 6._dp*g2p2*YuadjYuYu - 9._dp*TrYb3adjYb3*YuadjYuYu/10._dp -& 
&  9._dp*TrYuadjYu*YuadjYuYu - 9._dp*TrYw3adjYw3*YuadjYuYu/2._dp - 9._dp*TrYx3adjYx3*YuadjYuYu  
betaYu2 =  betaYu2- 4._dp*YuadjYuYuadjYuYu + 16._dp*g3p4*Yu*NGHg3 + 6._dp*g2p4*Yu*NGHw3 +& 
&  13._dp*g1p4*Yu*NGHx3/6._dp + 9._dp*g2p4*Yu*NGHx3/2._dp +& 
&  16._dp*g3p4*Yu*NGHx3/3._dp + 13._dp*g1p4*Yu*NGHxb3/6._dp +& 
&  9._dp*g2p4*Yu*NGHxb3/2._dp + 16._dp*g3p4*Yu*NGHxb3/3._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*TpYx3CYx3Yd - 7._dp*g1p2*Yd/15._dp - 3._dp*g2p2*Yd - 16._dp*g3p2*Yd/3._dp +& 
&  3._dp*TrYdadjYd*Yd + TrYeadjYe*Yd + 3._dp*YdadjYdYd + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = -2._dp*TpYx3CYx3TpYx3CYx3Yd + 2._dp*g1p2*TpYx3CYx3Yd + 6._dp*g2p2*TpYx3CYx3Yd -       & 
&  3._dp*TpYx3CYx3Yd*TrYb3adjYb3/5._dp - 6._dp*TpYx3CYx3Yd*TrYuadjYu - 3._dp*TpYx3CYx3Yd*TrYw3adjYw3 -& 
&  6._dp*TpYx3CYx3Yd*TrYx3adjYx3 + 287._dp*g1p4*Yd/90._dp + g1p2*g2p2*Yd +               & 
&  15._dp*g2p4*Yd/2._dp + 8._dp*g1p2*g3p2*Yd/9._dp + 8._dp*g2p2*g3p2*Yd - 16._dp*g3p4*Yd/9._dp -& 
&  6._dp*TrTpYx3CYx3YdadjYd*Yd - 2._dp*g1p2*TrYdadjYd*Yd/5._dp + 16._dp*g3p2*TrYdadjYd*Yd -& 
&  9._dp*TrYdadjYdYdadjYd*Yd - 3._dp*TrYdadjYuYuadjYd*Yd - 3._dp*TrYeadjYb3Yb3adjYe*Yd/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*Yd/5._dp - 3._dp*TrYeadjYeYeadjYe*Yd - 3._dp*TrYeadjYw3Yw3adjYe*Yd/2._dp -& 
&  2._dp*YdadjYdTpYx3CYx3Yd + 4._dp*g1p2*YdadjYdYd/5._dp + 6._dp*g2p2*YdadjYdYd -        & 
&  9._dp*TrYdadjYd*YdadjYdYd - 3._dp*TrYeadjYe*YdadjYdYd - 4._dp*YdadjYdYdadjYdYd +      & 
&  4._dp*g1p2*YdadjYuYu/5._dp - 3._dp*TrYb3adjYb3*YdadjYuYu/10._dp - 3._dp*TrYuadjYu*YdadjYuYu -& 
&  3._dp*TrYw3adjYw3*YdadjYuYu/2._dp - 3._dp*TrYx3adjYx3*YdadjYuYu - 2._dp*YdadjYuYuadjYdYd  
betaYd2 =  betaYd2- 2._dp*YdadjYuYuadjYuYu + 16._dp*g3p4*Yd*NGHg3 + 6._dp*g2p4*Yd*NGHw3 +& 
&  7._dp*g1p4*Yd*NGHx3/6._dp + 9._dp*g2p4*Yd*NGHx3/2._dp + & 
&  16._dp*g3p4*Yd*NGHx3/3._dp + 7._dp*g1p4*Yd*NGHxb3/6._dp +& 
&  9._dp*g2p4*Yd*NGHxb3/2._dp + 16._dp*g3p4*Yd*NGHxb3/3._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = -9._dp*g1p2*Ye/5._dp - 3._dp*g2p2*Ye + 3._dp*TrYdadjYd*Ye + TrYeadjYe*Ye + & 
&  3._dp*YeadjYb3Yb3/10._dp + 3._dp*YeadjYeYe + 3._dp*YeadjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = 27._dp*g1p4*Ye/2._dp + 9._dp*g1p2*g2p2*Ye/5._dp + 15._dp*g2p4*Ye/2._dp -              & 
&  6._dp*TrTpYx3CYx3YdadjYd*Ye - 2._dp*g1p2*TrYdadjYd*Ye/5._dp + 16._dp*g3p2*TrYdadjYd*Ye -& 
&  9._dp*TrYdadjYdYdadjYd*Ye - 3._dp*TrYdadjYuYuadjYd*Ye - 3._dp*TrYeadjYb3Yb3adjYe*Ye/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*Ye/5._dp - 3._dp*TrYeadjYeYeadjYe*Ye - 3._dp*TrYeadjYw3Yw3adjYe*Ye/2._dp -& 
&  9._dp*TrYb3adjYb3*YeadjYb3Yb3/100._dp - 9._dp*TrYuadjYu*YeadjYb3Yb3/10._dp -          & 
&  9._dp*TrYw3adjYw3*YeadjYb3Yb3/20._dp - 9._dp*TrYx3adjYx3*YeadjYb3Yb3/10._dp -         & 
&  9._dp*YeadjYb3Yb3adjYb3Yb3/50._dp - 3._dp*YeadjYb3Yb3adjYeYe/5._dp - 3._dp*YeadjYb3Yb3adjYw3Yw3/10._dp +& 
&  6._dp*g2p2*YeadjYeYe - 9._dp*TrYdadjYd*YeadjYeYe - 3._dp*TrYeadjYe*YeadjYeYe -        & 
&  4._dp*YeadjYeYeadjYeYe + 6._dp*g2p2*YeadjYw3Yw3 - 9._dp*TrYb3adjYb3*YeadjYw3Yw3/20._dp -& 
&  9._dp*TrYuadjYu*YeadjYw3Yw3/2._dp - 9._dp*TrYw3adjYw3*YeadjYw3Yw3/4._dp -             & 
&  9._dp*TrYx3adjYx3*YeadjYw3Yw3/2._dp - 9._dp*YeadjYw3Yw3adjYb3Yb3/40._dp  
betaYe2 =  betaYe2- 3._dp*YeadjYw3Yw3adjYeYe - 3._dp*YeadjYw3Yw3adjYw3Yw3/2._dp + 6._dp*g2p4*Ye*NGHw3 +& 
&  9._dp*g1p4*Ye*NGHx3/2._dp + 9._dp*g2p4*Ye*NGHx3/2._dp + & 
&  9._dp*g1p4*Ye*NGHxb3/2._dp + 9._dp*g2p4*Ye*NGHxb3/2._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = -3._dp*g1p2*Yb3/5._dp - 3._dp*g2p2*Yb3 + 3._dp*TrYb3adjYb3*Yb3/10._dp +   & 
&  3._dp*TrYuadjYu*Yb3 + 3._dp*TrYw3adjYw3*Yb3/2._dp + 3._dp*TrYx3adjYx3*Yb3 +           & 
&  9._dp*Yb3adjYb3Yb3/10._dp + Yb3adjYeYe + 3._dp*Yb3adjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = 207._dp*g1p4*Yb3/50._dp + 9._dp*g1p2*g2p2*Yb3/5._dp + 15._dp*g2p4*Yb3/2._dp -         & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yb3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yb3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yb3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yb3 - 3._dp*TrYuadjYdYdadjYu*Yb3 +& 
&  4._dp*g1p2*TrYuadjYu*Yb3/5._dp + 16._dp*g3p2*TrYuadjYu*Yb3 - 9._dp*TrYuadjYuYuadjYu*Yb3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yb3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yb3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yb3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yb3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yb3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yb3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yb3 + 9._dp*g1p2*Yb3adjYb3Yb3/25._dp +& 
&  9._dp*g2p2*Yb3adjYb3Yb3/5._dp - 27._dp*TrYb3adjYb3*Yb3adjYb3Yb3/100._dp -             & 
&  27._dp*TrYuadjYu*Yb3adjYb3Yb3/10._dp - 27._dp*TrYw3adjYw3*Yb3adjYb3Yb3/20._dp -       & 
&  27._dp*TrYx3adjYx3*Yb3adjYb3Yb3/10._dp - 9._dp*Yb3adjYb3Yb3adjYb3Yb3/25._dp -         & 
&  3._dp*Yb3adjYb3Yb3adjYw3Yw3/10._dp + 6._dp*g1p2*Yb3adjYeYe/5._dp - 3._dp*TrYdadjYd*Yb3adjYeYe  
betaYb32 =  betaYb32- TrYeadjYe*Yb3adjYeYe - 3._dp*Yb3adjYeYeadjYb3Yb3/5._dp - 2._dp*Yb3adjYeYeadjYeYe +    & 
&  6._dp*g2p2*Yb3adjYw3Yw3 - 9._dp*TrYb3adjYb3*Yb3adjYw3Yw3/20._dp - 9._dp*TrYuadjYu*Yb3adjYw3Yw3/2._dp -& 
&  9._dp*TrYw3adjYw3*Yb3adjYw3Yw3/4._dp - 9._dp*TrYx3adjYx3*Yb3adjYw3Yw3/2._dp -         & 
&  9._dp*Yb3adjYw3Yw3adjYb3Yb3/8._dp - 3._dp*Yb3adjYw3Yw3adjYw3Yw3/2._dp +               & 
&  6._dp*g2p4*Yb3*NGHw3 + 3._dp*g1p4*Yb3*NGHx3/2._dp +     & 
&  9._dp*g2p4*Yb3*NGHx3/2._dp + 3._dp*g1p4*Yb3*NGHxb3/2._dp +& 
&  9._dp*g2p4*Yb3*NGHxb3/2._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = -3._dp*g1p2*Yw3/5._dp - 7._dp*g2p2*Yw3 + 3._dp*TrYb3adjYb3*Yw3/10._dp +   & 
&  3._dp*TrYuadjYu*Yw3 + 3._dp*TrYw3adjYw3*Yw3/2._dp + 3._dp*TrYx3adjYx3*Yw3 +           & 
&  3._dp*Yw3adjYb3Yb3/10._dp + Yw3adjYeYe + 5._dp*Yw3adjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = 207._dp*g1p4*Yw3/50._dp + 9._dp*g1p2*g2p2*Yw3/5._dp + 55._dp*g2p4*Yw3/2._dp -         & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yw3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yw3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yw3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yw3 - 3._dp*TrYuadjYdYdadjYu*Yw3 +& 
&  4._dp*g1p2*TrYuadjYu*Yw3/5._dp + 16._dp*g3p2*TrYuadjYu*Yw3 - 9._dp*TrYuadjYuYuadjYu*Yw3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yw3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yw3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yw3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yw3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yw3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yw3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yw3 - 9._dp*TrYb3adjYb3*Yw3adjYb3Yb3/100._dp -& 
&  9._dp*TrYuadjYu*Yw3adjYb3Yb3/10._dp - 9._dp*TrYw3adjYw3*Yw3adjYb3Yb3/20._dp -         & 
&  9._dp*TrYx3adjYx3*Yw3adjYb3Yb3/10._dp - 9._dp*Yw3adjYb3Yb3adjYb3Yb3/50._dp -          & 
&  3._dp*Yw3adjYb3Yb3adjYw3Yw3/5._dp + 6._dp*g1p2*Yw3adjYeYe/5._dp - 3._dp*TrYdadjYd*Yw3adjYeYe -& 
&  TrYeadjYe*Yw3adjYeYe - 2._dp*Yw3adjYeYeadjYeYe - Yw3adjYeYeadjYw3Yw3 + 3._dp*g1p2*Yw3adjYw3Yw3/5._dp  
betaYw32 =  betaYw32+ 5._dp*g2p2*Yw3adjYw3Yw3 - 3._dp*TrYb3adjYb3*Yw3adjYw3Yw3/4._dp - 15._dp*TrYuadjYu*Yw3adjYw3Yw3/2._dp -& 
&  15._dp*TrYw3adjYw3*Yw3adjYw3Yw3/4._dp - 15._dp*TrYx3adjYx3*Yw3adjYw3Yw3/2._dp -       & 
&  9._dp*Yw3adjYw3Yw3adjYb3Yb3/40._dp - 3._dp*Yw3adjYw3Yw3adjYw3Yw3 + 14._dp*g2p4*Yw3*NGHw3 +& 
&  3._dp*g1p4*Yw3*NGHx3/2._dp + 21._dp*g2p4*Yw3*NGHx3/2._dp +& 
&  3._dp*g1p4*Yw3*NGHxb3/2._dp + 21._dp*g2p4*Yw3*NGHxb3/2._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = -19._dp*g1p2*Yx3/15._dp - 3._dp*g2p2*Yx3 - 16._dp*g3p2*Yx3/3._dp +        & 
&  3._dp*TrYb3adjYb3*Yx3/10._dp + 3._dp*TrYuadjYu*Yx3 + 3._dp*TrYw3adjYw3*Yx3/2._dp +    & 
&  3._dp*TrYx3adjYx3*Yx3 + 3._dp*Yx3adjYx3Yx3 + 2._dp*Yx3CYdTpYd

 
 
If (TwoLoopRGE) Then 
betaYx32 = 4123._dp*g1p4*Yx3/450._dp + 17._dp*g1p2*g2p2*Yx3/5._dp + 15._dp*g2p4*Yx3/2._dp +      & 
&  232._dp*g1p2*g3p2*Yx3/45._dp + 8._dp*g2p2*g3p2*Yx3 - 16._dp*g3p4*Yx3/9._dp -          & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yx3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yx3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yx3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yx3 - 3._dp*TrYuadjYdYdadjYu*Yx3 +& 
&  4._dp*g1p2*TrYuadjYu*Yx3/5._dp + 16._dp*g3p2*TrYuadjYu*Yx3 - 9._dp*TrYuadjYuYuadjYu*Yx3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yx3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yx3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yx3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yx3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yx3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yx3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yx3 + 8._dp*g1p2*Yx3adjYx3Yx3/5._dp +& 
&  6._dp*g2p2*Yx3adjYx3Yx3 - 9._dp*TrYb3adjYb3*Yx3adjYx3Yx3/10._dp - 9._dp*TrYuadjYu*Yx3adjYx3Yx3 -& 
&  9._dp*TrYw3adjYw3*Yx3adjYx3Yx3/2._dp - 9._dp*TrYx3adjYx3*Yx3adjYx3Yx3 -               & 
&  4._dp*Yx3adjYx3Yx3adjYx3Yx3 + 2._dp*g1p2*Yx3CYdTpYd/5._dp + 6._dp*g2p2*Yx3CYdTpYd  
betaYx32 =  betaYx32- 6._dp*TrYdadjYd*Yx3CYdTpYd - 2._dp*TrYeadjYe*Yx3CYdTpYd - 2._dp*Yx3CYdTpYdadjYx3Yx3 - & 
&  2._dp*Yx3CYdTpYdCYdTpYd - 2._dp*Yx3CYdTpYuCYuTpYd + 16._dp*g3p4*Yx3*NGHg3 +& 
&  6._dp*g2p4*Yx3*NGHw3 + 19._dp*g1p4*Yx3*NGHx3/6._dp +    & 
&  9._dp*g2p4*Yx3*NGHx3/2._dp + 16._dp*g3p4*Yx3*NGHx3/3._dp +& 
&  19._dp*g1p4*Yx3*NGHxb3/6._dp + 9._dp*g2p4*Yx3*NGHxb3/2._dp +& 
&  16._dp*g3p4*Yx3*NGHxb3/3._dp

 
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
 
  Integer i1, i2, i3, i4, SumI 
 
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

  Integer i1, i2, i3, i4, SumI 
 
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

Subroutine rge555(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2,i3,i4 
Integer :: j1,j2,j3,j4,j5,j6,j7 
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
& ,MassB,betaMassB1,betaMassB2,DMassB,MassWB,betaMassWB1,betaMassWB2,DMassWB,            & 
& MassG,betaMassG1,betaMassG2,DMassG
Complex(dp) :: Tr1(3),Tr2(3),Tr3(3) 
Complex(dp) :: md2adjYx3(3,3),md2CYd(3,3),me2CYe(3,3),mHb32CYb3(3,3),mHw32CYw3(3,3),mHxb32CYx3(3,3)  & 
& ,ml2adjYb3(3,3),ml2adjYe(3,3),ml2adjYw3(3,3),mq2adjYd(3,3),mq2adjYu(3,3)               & 
& ,mu2CYu(3,3),Yb3adjYb3(3,3),YdadjYd(3,3),YeadjYe(3,3),YuadjYu(3,3),Yw3adjYw3(3,3)      & 
& ,Yx3adjYx3(3,3),adjYb3MBM3(3,3),adjYb3mHb32(3,3),adjYb3Yb3(3,3),adjYb3BMBM3(3,3)       & 
& ,adjYb3TYb3(3,3),adjYdmd2(3,3),adjYdYd(3,3),adjYdTYd(3,3),adjYeme2(3,3),               & 
& adjYeYe(3,3),adjYeTYe(3,3),adjYumu2(3,3),adjYuYu(3,3),adjYuTYu(3,3),adjYw3mHw32(3,3)   & 
& ,adjYw3MWM3(3,3),adjYw3Yw3(3,3),adjYw3BMWM3(3,3),adjYw3TYw3(3,3),adjYx3mHxb32(3,3)     & 
& ,adjYx3Yx3(3,3),adjYx3TYx3(3,3),CYb3ml2(3,3),CYb3TpYb3(3,3),CYb3TpTYb3(3,3)            & 
& ,CYdmq2(3,3),CYdTpYd(3,3),CYdTpTYd(3,3),CYeml2(3,3),CYumq2(3,3),CYw3ml2(3,3)           & 
& ,CYw3TpYw3(3,3),CYw3TpTYw3(3,3),CYx3md2(3,3),CYx3Yd(3,3),CYx3TYd(3,3),CYx3TpYx3(3,3)   & 
& ,CYx3TpTYx3(3,3),TYb3adjYb3(3,3),TYb3adjTYb3(3,3),TYdadjYd(3,3),TYdadjTYd(3,3)         & 
& ,TYeadjYe(3,3),TYeadjTYe(3,3),TYuadjYu(3,3),TYuadjTYu(3,3),TYw3adjYw3(3,3)             & 
& ,TYw3adjTYw3(3,3),TYx3adjYx3(3,3),TYx3adjTYx3(3,3),TpYb3CYb3(3,3),TpYdCYd(3,3)         & 
& ,TpYeCYe(3,3),TpYuCYu(3,3),TpYw3CYw3(3,3),TpYx3CYx3(3,3),TpTYb3CTYb3(3,3)              & 
& ,TpTYdCTYd(3,3),TpTYeCTYe(3,3),TpTYuCTYu(3,3),TpTYw3CTYw3(3,3),TpTYx3CTYx3(3,3)        & 
& ,MBM3CYb3TpYb3(3,3),MBM3CYb3TpTYb3(3,3),md2YdadjYd(3,3),md2TpYx3CYx3(3,3)              & 
& ,me2YeadjYe(3,3),mHb32Yb3adjYb3(3,3),mHw32Yw3adjYw3(3,3),mHxb32Yx3adjYx3(3,3)          & 
& ,ml2TpYb3CYb3(3,3),ml2TpYeCYe(3,3),ml2TpYw3CYw3(3,3),mq2TpYdCYd(3,3),mq2TpYuCYu(3,3)   & 
& ,mu2YuadjYu(3,3),MWM3CYw3TpYw3(3,3),MWM3CYw3TpTYw3(3,3),MXM3CYx3TpYx3(3,3)             & 
& ,MXM3CYx3TpTYx3(3,3),Yb3ml2adjYb3(3,3),Yb3adjYb3MBM3(3,3),Yb3adjYb3mHb32(3,3)          & 
& ,Yb3adjYb3Yb3(3,3),Yb3adjYb3BMBM3(3,3),Yb3adjYb3TYb3(3,3),Yb3adjYeYe(3,3)              & 
& ,Yb3adjYeTYe(3,3),Yb3adjYw3Yw3(3,3),Yb3adjYw3TYw3(3,3),Ydmq2adjYd(3,3),YdadjYdmd2(3,3) & 
& ,YdadjYdYd(3,3),YdadjYdTYd(3,3),YdadjYuYu(3,3),YdadjYuTYu(3,3),Yeml2adjYe(3,3)         & 
& ,YeadjYb3Yb3(3,3),YeadjYb3TYb3(3,3),YeadjYeme2(3,3),YeadjYeYe(3,3),YeadjYeTYe(3,3)     & 
& ,YeadjYw3Yw3(3,3),YeadjYw3TYw3(3,3),Yumq2adjYu(3,3),YuadjYdYd(3,3),YuadjYdTYd(3,3)     & 
& ,YuadjYumu2(3,3),YuadjYuYu(3,3),YuadjYuTYu(3,3),Yw3ml2adjYw3(3,3),Yw3adjYb3Yb3(3,3)    & 
& ,Yw3adjYb3TYb3(3,3),Yw3adjYeYe(3,3),Yw3adjYeTYe(3,3),Yw3adjYw3mHw32(3,3)               & 
& ,Yw3adjYw3MWM3(3,3),Yw3adjYw3Yw3(3,3),Yw3adjYw3BMWM3(3,3),Yw3adjYw3TYw3(3,3)           & 
& ,Yx3md2adjYx3(3,3),Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3(3,3),Yx3adjYx3TYx3(3,3)           & 
& ,Yx3CYdTpYd(3,3),Yx3CYdTpTYd(3,3),BMBM3CYb3TpYb3(3,3),BMWM3CYw3TpYw3(3,3)              & 
& ,BMXM3CYx3TpYx3(3,3),TYb3adjYb3MBM3(3,3),TYb3adjYb3Yb3(3,3),TYb3adjYeYe(3,3)           & 
& ,TYb3adjYw3Yw3(3,3),TYdadjYdYd(3,3),TYdadjYuYu(3,3),TYeadjYb3Yb3(3,3),TYeadjYeYe(3,3)  & 
& ,TYeadjYw3Yw3(3,3),TYuadjYdYd(3,3),TYuadjYuYu(3,3),TYw3adjYb3Yb3(3,3),TYw3adjYeYe(3,3) & 
& ,TYw3adjYw3MWM3(3,3),TYw3adjYw3Yw3(3,3),TYx3adjYx3Yx3(3,3),TYx3CYdTpYd(3,3)            & 
& ,TpYb3mHb32CYb3(3,3),TpYb3CYb3ml2(3,3),TpYdmd2CYd(3,3),TpYdCYdmq2(3,3),TpYeme2CYe(3,3) & 
& ,TpYeCYeml2(3,3),TpYumu2CYu(3,3),TpYuCYumq2(3,3),TpYw3mHw32CYw3(3,3),TpYw3CYw3ml2(3,3) & 
& ,TpYx3mHxb32CYx3(3,3),TpYx3CYx3md2(3,3),TpYx3CYx3Yd(3,3),TpYx3CYx3TYd(3,3)             

Complex(dp) :: TpTYx3CYx3Yd(3,3)

Complex(dp) :: mHxb32Yx3(3,3),Yb3adjYe(3,3),Yb3adjYw3(3,3),Yb3adjTYb3(3,3),Yb3adjTYe(3,3)            & 
& ,Yb3adjTYw3(3,3),YdadjYu(3,3),YdadjTYd(3,3),YdadjTYu(3,3),YeadjYb3(3,3),               & 
& YeadjYw3(3,3),YeadjTYb3(3,3),YeadjTYe(3,3),YeadjTYw3(3,3),YuadjYd(3,3),YuadjTYd(3,3)   & 
& ,YuadjTYu(3,3),Yw3adjYb3(3,3),Yw3adjYe(3,3),Yw3adjTYb3(3,3),Yw3adjTYe(3,3)             & 
& ,Yw3adjTYw3(3,3),Yx3adjTYx3(3,3),Yx3CYd(3,3),Yx3CTYd(3,3),adjYdadjTYx3(3,3)            & 
& ,adjYdTpYx3(3,3),adjYdTpTYx3(3,3),adjTYdadjYx3(3,3),CYb3TpYw3(3,3),CYb3TpTYw3(3,3)     & 
& ,CYeTpYb3(3,3),CYeTpYw3(3,3),CYeTpTYb3(3,3),CYeTpTYw3(3,3),CYuTpYd(3,3),               & 
& CYuTpTYd(3,3),CYw3TpYb3(3,3),CYw3TpTYb3(3,3),CTYb3adjYb3(3,3),CTYb3adjYe(3,3)          & 
& ,CTYb3adjYw3(3,3),CTYdadjYd(3,3),CTYdadjYu(3,3),CTYeadjYb3(3,3),CTYeadjYe(3,3)         & 
& ,CTYeadjYw3(3,3),CTYuadjYd(3,3),CTYuadjYu(3,3),CTYw3adjYb3(3,3),CTYw3adjYe(3,3)        & 
& ,CTYw3adjYw3(3,3),CTYx3adjYx3(3,3),TYb3adjYe(3,3),TYb3adjYw3(3,3),TYb3adjTYe(3,3)      & 
& ,TYb3adjTYw3(3,3),TYdadjYu(3,3),TYdadjTYu(3,3),TYeadjYb3(3,3),TYeadjYw3(3,3)           & 
& ,TYeadjTYb3(3,3),TYeadjTYw3(3,3),TYuadjYd(3,3),TYuadjTYd(3,3),TYw3adjYb3(3,3)          & 
& ,TYw3adjYe(3,3),TYw3adjTYb3(3,3),TYw3adjTYe(3,3),TYx3CTYd(3,3),TpYb3CTYb3(3,3)         & 
& ,TpYdadjYx3(3,3),TpYdadjTYx3(3,3),TpYdCTYd(3,3),TpYeCTYe(3,3),TpYuCTYu(3,3)            & 
& ,TpYw3CTYw3(3,3),TpYx3CTYx3(3,3),TpTYb3CYb3(3,3),TpTYdadjYx3(3,3),TpTYdCYd(3,3)        & 
& ,TpTYeCYe(3,3),TpTYuCYu(3,3),TpTYw3CYw3(3,3),TpTYx3CYx3(3,3),md2YdadjYu(3,3)           & 
& ,me2YeadjYb3(3,3),me2YeadjYw3(3,3),mHb32Yb3adjYe(3,3),mHb32Yb3adjYw3(3,3)              & 
& ,mHw32Yw3adjYb3(3,3),mHw32Yw3adjYe(3,3),mHxb32Yx3CYd(3,3),mq2TpYdadjYx3(3,3)           & 
& ,mu2YuadjYd(3,3),Yb3ml2adjYe(3,3),Yb3ml2adjYw3(3,3),Yb3adjYeme2(3,3),Yb3adjYw3mHw32(3,3)& 
& ,Yb3adjYw3MWM3(3,3),Yb3adjYw3BMWM3(3,3),Ydmq2adjYu(3,3),YdadjYdTpYx3(3,3)              & 
& ,YdadjYdTpTYx3(3,3),YdadjYumu2(3,3),YdadjTYdadjYx3(3,3),Yeml2adjYb3(3,3)               & 
& ,Yeml2adjYw3(3,3),YeadjYb3MBM3(3,3),YeadjYb3mHb32(3,3),YeadjYb3BMBM3(3,3)              & 
& ,YeadjYw3mHw32(3,3),YeadjYw3MWM3(3,3),YeadjYw3BMWM3(3,3),Yumq2adjYd(3,3)               & 
& ,YuadjYdmd2(3,3),Yw3ml2adjYb3(3,3),Yw3ml2adjYe(3,3),Yw3adjYb3MBM3(3,3),Yw3adjYb3mHb32(3,3)& 
& ,Yw3adjYb3BMBM3(3,3),Yw3adjYeme2(3,3),Yx3md2CYd(3,3),Yx3CYdmq2(3,3),adjYb3Yb3adjYb3(3,3)& 
& ,adjYb3Yb3adjYe(3,3),adjYb3Yb3adjYw3(3,3),adjYb3Yb3adjTYb3(3,3),adjYb3Yb3adjTYe(3,3)   & 
& ,adjYb3Yb3adjTYw3(3,3),adjYb3TYb3adjYb3(3,3),adjYb3TYb3adjYe(3,3),adjYb3TYb3adjYw3(3,3)& 
& ,adjYb3TYb3adjTYb3(3,3),adjYb3TYb3adjTYe(3,3),adjYb3TYb3adjTYw3(3,3),adjYdYdadjYd(3,3) & 
& ,adjYdYdadjYu(3,3),adjYdYdadjTYd(3,3),adjYdYdadjTYu(3,3),adjYdadjYx3Yx3(3,3)           & 
& ,adjYdTYdadjYd(3,3),adjYdTYdadjYu(3,3),adjYdTYdadjTYd(3,3),adjYdTYdadjTYu(3,3)         & 
& ,adjYdTpYx3CYx3(3,3),adjYdTpYx3CTYx3(3,3),adjYdTpTYx3CTYx3(3,3),adjYeYeadjYb3(3,3)     & 
& ,adjYeYeadjYe(3,3),adjYeYeadjYw3(3,3),adjYeYeadjTYb3(3,3),adjYeYeadjTYe(3,3)           & 
& ,adjYeYeadjTYw3(3,3),adjYeTYeadjYb3(3,3),adjYeTYeadjYe(3,3),adjYeTYeadjYw3(3,3)        & 
& ,adjYeTYeadjTYb3(3,3),adjYeTYeadjTYe(3,3),adjYeTYeadjTYw3(3,3),adjYuYuadjYd(3,3)       & 
& ,adjYuYuadjYu(3,3),adjYuYuadjTYd(3,3),adjYuYuadjTYu(3,3),adjYuTYuadjYd(3,3)            & 
& ,adjYuTYuadjYu(3,3),adjYuTYuadjTYd(3,3),adjYuTYuadjTYu(3,3),adjYw3Yw3adjYb3(3,3)       

Complex(dp) :: adjYw3Yw3adjYe(3,3),adjYw3Yw3adjYw3(3,3),adjYw3Yw3adjTYb3(3,3),adjYw3Yw3adjTYe(3,3)   & 
& ,adjYw3Yw3adjTYw3(3,3),adjYw3TYw3adjYb3(3,3),adjYw3TYw3adjYe(3,3),adjYw3TYw3adjYw3(3,3)& 
& ,adjYw3TYw3adjTYb3(3,3),adjYw3TYw3adjTYe(3,3),adjYw3TYw3adjTYw3(3,3),adjYx3mHxb32Yx3(3,3)& 
& ,adjYx3Yx3adjYx3(3,3),adjYx3Yx3adjTYx3(3,3),adjYx3Yx3CYd(3,3),adjYx3Yx3CTYd(3,3)       & 
& ,adjYx3TYx3adjYx3(3,3),adjYx3TYx3adjTYx3(3,3),adjYx3TYx3CTYd(3,3),adjTYb3TYb3adjYb3(3,3)& 
& ,adjTYb3TYb3adjYe(3,3),adjTYb3TYb3adjYw3(3,3),adjTYdTYdadjYd(3,3),adjTYdTYdadjYu(3,3)  & 
& ,adjTYeTYeadjYb3(3,3),adjTYeTYeadjYe(3,3),adjTYeTYeadjYw3(3,3),adjTYuTYuadjYd(3,3)     & 
& ,adjTYuTYuadjYu(3,3),adjTYw3TYw3adjYb3(3,3),adjTYw3TYw3adjYe(3,3),adjTYw3TYw3adjYw3(3,3)& 
& ,adjTYx3TYx3adjYx3(3,3),CYb3TpYb3CYb3(3,3),CYb3TpYb3CTYb3(3,3),CYb3TpYw3CYw3(3,3)      & 
& ,CYb3TpYw3CTYw3(3,3),CYb3TpTYb3CTYb3(3,3),CYb3TpTYw3CTYw3(3,3),CYdTpYdadjYx3(3,3)      & 
& ,CYdTpYdadjTYx3(3,3),CYdTpYdCYd(3,3),CYdTpYdCTYd(3,3),CYdTpTYdCTYd(3,3),               & 
& CYeTpYeCYe(3,3),CYeTpYeCTYe(3,3),CYeTpTYeCTYe(3,3),CYuTpYuCYu(3,3),CYuTpYuCTYu(3,3)    & 
& ,CYuTpTYuCTYu(3,3),CYw3TpYb3CYb3(3,3),CYw3TpYb3CTYb3(3,3),CYw3TpYw3CYw3(3,3)           & 
& ,CYw3TpYw3CTYw3(3,3),CYw3TpTYb3CTYb3(3,3),CYw3TpTYw3CTYw3(3,3),CYx3YdadjYd(3,3)        & 
& ,CYx3TYdadjYd(3,3),CYx3TpYx3CYx3(3,3),CYx3TpYx3CTYx3(3,3),CYx3TpTYx3CTYx3(3,3)         & 
& ,CTYb3TpTYb3CYb3(3,3),CTYb3TpTYw3CYw3(3,3),CTYdTpTYdadjYx3(3,3),CTYdTpTYdCYd(3,3)      & 
& ,CTYeTpTYeCYe(3,3),CTYuTpTYuCYu(3,3),CTYw3TpTYb3CYb3(3,3),CTYw3TpTYw3CYw3(3,3)         & 
& ,CTYx3TpTYx3CYx3(3,3),TYb3adjYw3MWM3(3,3),TYb3TpYb3CYb3(3,3),TYb3TpYw3CYw3(3,3)        & 
& ,TYdadjYdadjTYx3(3,3),TYdadjYdTpYx3(3,3),TYdTpYdCYd(3,3),TYeadjYb3MBM3(3,3)            & 
& ,TYeadjYw3MWM3(3,3),TYeTpYeCYe(3,3),TYuTpYuCYu(3,3),TYw3adjYb3MBM3(3,3),               & 
& TYw3TpYb3CYb3(3,3),TYw3TpYw3CYw3(3,3),TYx3TpYx3CYx3(3,3),TpYb3CYb3TpYb3(3,3)           & 
& ,TpYb3CYb3TpYw3(3,3),TpYb3CYb3TpTYb3(3,3),TpYb3CYb3TpTYw3(3,3),TpYb3CTYb3adjYb3(3,3)   & 
& ,TpYb3CTYb3adjYe(3,3),TpYb3CTYb3adjYw3(3,3),TpYdmd2adjYx3(3,3),TpYdadjYx3mHxb32(3,3)   & 
& ,TpYdadjYx3Yx3(3,3),TpYdadjYx3TYx3(3,3),TpYdCYdTpYd(3,3),TpYdCYdTpTYd(3,3)             & 
& ,TpYdCTYdadjYd(3,3),TpYdCTYdadjYu(3,3),TpYeCYeTpYb3(3,3),TpYeCYeTpYw3(3,3)             & 
& ,TpYeCYeTpTYb3(3,3),TpYeCYeTpTYw3(3,3),TpYeCTYeadjYb3(3,3),TpYeCTYeadjYe(3,3)          & 
& ,TpYeCTYeadjYw3(3,3),TpYuCYuTpYd(3,3),TpYuCYuTpTYd(3,3),TpYuCTYuadjYd(3,3)             & 
& ,TpYuCTYuadjYu(3,3),TpYw3CYw3TpYb3(3,3),TpYw3CYw3TpYw3(3,3),TpYw3CYw3TpTYb3(3,3)       & 
& ,TpYw3CYw3TpTYw3(3,3),TpYw3CTYw3adjYb3(3,3),TpYw3CTYw3adjYe(3,3),TpYw3CTYw3adjYw3(3,3) & 
& ,TpYx3CYx3TpYx3(3,3),TpYx3CYx3TpTYx3(3,3),TpYx3CTYx3adjYx3(3,3),TpTYb3CYb3TpYb3(3,3)   & 
& ,TpTYb3CYb3TpYw3(3,3),TpTYb3CTYb3adjYb3(3,3),TpTYb3CTYb3adjYe(3,3),TpTYb3CTYb3adjYw3(3,3)& 
& ,TpTYdCYdTpYd(3,3),TpTYdCTYdadjYd(3,3),TpTYdCTYdadjYu(3,3),TpTYeCYeTpYb3(3,3)          & 
& ,TpTYeCYeTpYw3(3,3),TpTYeCTYeadjYb3(3,3),TpTYeCTYeadjYe(3,3),TpTYeCTYeadjYw3(3,3)      & 
& ,TpTYuCYuTpYd(3,3),TpTYuCTYuadjYd(3,3),TpTYuCTYuadjYu(3,3),TpTYw3CYw3TpYb3(3,3)        & 
& ,TpTYw3CYw3TpYw3(3,3),TpTYw3CTYw3adjYb3(3,3),TpTYw3CTYw3adjYe(3,3),TpTYw3CTYw3adjYw3(3,3)& 
& ,TpTYx3CYx3TpYx3(3,3),TpTYx3CTYx3adjYx3(3,3),md2adjYx3Yx3adjYx3(3,3),md2adjYx3Yx3CYd(3,3)& 
& ,md2CYdTpYdadjYx3(3,3),md2CYdTpYdCYd(3,3),me2CYeTpYeCYe(3,3),mHb32CYb3TpYb3CYb3(3,3)   

Complex(dp) :: mHb32CYb3TpYw3CYw3(3,3),mHw32CYw3TpYb3CYb3(3,3),mHw32CYw3TpYw3CYw3(3,3)               & 
& ,mHxb32CYx3YdadjYd(3,3),mHxb32CYx3TpYx3CYx3(3,3),ml2adjYb3Yb3adjYb3(3,3)               & 
& ,ml2adjYb3Yb3adjYe(3,3),ml2adjYb3Yb3adjYw3(3,3),ml2adjYeYeadjYb3(3,3),ml2adjYeYeadjYe(3,3)& 
& ,ml2adjYeYeadjYw3(3,3),ml2adjYw3Yw3adjYb3(3,3),ml2adjYw3Yw3adjYe(3,3),ml2adjYw3Yw3adjYw3(3,3)& 
& ,mq2adjYdYdadjYd(3,3),mq2adjYdYdadjYu(3,3),mq2adjYdTpYx3CYx3(3,3),mq2adjYuYuadjYd(3,3) & 
& ,mq2adjYuYuadjYu(3,3),mq2adjYx3Yx3CYd(3,3),mu2CYuTpYuCYu(3,3),Yb3adjYb3Yb3adjYb3(3,3)  & 
& ,Yb3adjYb3TYb3adjYb3(3,3),Yb3adjYb3TYb3adjTYb3(3,3),Yb3adjYeYeadjYb3(3,3)              & 
& ,Yb3adjYeTYeadjYb3(3,3),Yb3adjYeTYeadjTYb3(3,3),Yb3adjYw3Yw3adjYb3(3,3),               & 
& Yb3adjYw3TYw3adjYb3(3,3),Yb3adjYw3TYw3adjTYb3(3,3),Yb3adjTYb3TYb3adjYb3(3,3)           & 
& ,Yb3adjTYeTYeadjYb3(3,3),Yb3adjTYw3TYw3adjYb3(3,3),Yb3TpTYb3CTYb3adjYb3(3,3)           & 
& ,Yb3TpTYeCTYeadjYb3(3,3),Yb3TpTYw3CTYw3adjYb3(3,3),YdadjYdYdadjYd(3,3),YdadjYdTYdadjYd(3,3)& 
& ,YdadjYdTYdadjTYd(3,3),YdadjYdTpYx3CYx3(3,3),YdadjYdTpTYx3CTYx3(3,3),YdadjYuYuadjYd(3,3)& 
& ,YdadjYuTYuadjYd(3,3),YdadjYuTYuadjTYd(3,3),YdadjTYdTYdadjYd(3,3),YdadjTYuTYuadjYd(3,3)& 
& ,YdTpTYdCTYdadjYd(3,3),YdTpTYuCTYuadjYd(3,3),YeadjYb3Yb3adjYe(3,3),YeadjYb3TYb3adjYe(3,3)& 
& ,YeadjYb3TYb3adjTYe(3,3),YeadjYeYeadjYe(3,3),YeadjYeTYeadjYe(3,3),YeadjYeTYeadjTYe(3,3)& 
& ,YeadjYw3Yw3adjYe(3,3),YeadjYw3TYw3adjYe(3,3),YeadjYw3TYw3adjTYe(3,3),YeadjTYb3TYb3adjYe(3,3)& 
& ,YeadjTYeTYeadjYe(3,3),YeadjTYw3TYw3adjYe(3,3),YeTpTYb3CTYb3adjYe(3,3),YeTpTYeCTYeadjYe(3,3)& 
& ,YeTpTYw3CTYw3adjYe(3,3),YuadjYdYdadjYu(3,3),YuadjYdTYdadjYu(3,3),YuadjYdTYdadjTYu(3,3)& 
& ,YuadjYuYuadjYu(3,3),YuadjYuTYuadjYu(3,3),YuadjYuTYuadjTYu(3,3),YuadjTYdTYdadjYu(3,3)  & 
& ,YuadjTYuTYuadjYu(3,3),YuTpTYdCTYdadjYu(3,3),YuTpTYuCTYuadjYu(3,3),Yw3adjYb3Yb3adjYw3(3,3)& 
& ,Yw3adjYb3TYb3adjYw3(3,3),Yw3adjYb3TYb3adjTYw3(3,3),Yw3adjYeYeadjYw3(3,3)              & 
& ,Yw3adjYeTYeadjYw3(3,3),Yw3adjYeTYeadjTYw3(3,3),Yw3adjYw3Yw3adjYw3(3,3),               & 
& Yw3adjYw3TYw3adjYw3(3,3),Yw3adjYw3TYw3adjTYw3(3,3),Yw3adjTYb3TYb3adjYw3(3,3)           & 
& ,Yw3adjTYeTYeadjYw3(3,3),Yw3adjTYw3TYw3adjYw3(3,3),Yw3TpTYb3CTYb3adjYw3(3,3)           & 
& ,Yw3TpTYeCTYeadjYw3(3,3),Yw3TpTYw3CTYw3adjYw3(3,3),Yx3adjYx3Yx3adjYx3(3,3)             & 
& ,Yx3adjYx3TYx3adjYx3(3,3),Yx3adjYx3TYx3adjTYx3(3,3),Yx3adjTYx3TYx3adjYx3(3,3)          & 
& ,Yx3CYdTpYdadjYx3(3,3),Yx3CTYdTpTYdadjYx3(3,3),Yx3TYdadjYdadjTYx3(3,3),Yx3TpTYx3CTYx3adjYx3(3,3)& 
& ,adjYb3mHb32Yb3adjYb3(3,3),adjYb3mHb32Yb3adjYe(3,3),adjYb3mHb32Yb3adjYw3(3,3)          & 
& ,adjYb3Yb3ml2adjYb3(3,3),adjYb3Yb3ml2adjYe(3,3),adjYb3Yb3ml2adjYw3(3,3),               & 
& adjYb3Yb3adjYb3MBM3(3,3),adjYb3Yb3adjYb3mHb32(3,3),adjYb3Yb3adjYb3Yb3(3,3)             & 
& ,adjYb3Yb3adjYb3BMBM3(3,3),adjYb3Yb3adjYb3TYb3(3,3),adjYb3Yb3adjYeme2(3,3)             & 
& ,adjYb3Yb3adjYeYe(3,3),adjYb3Yb3adjYeTYe(3,3),adjYb3Yb3adjYw3mHw32(3,3),               & 
& adjYb3Yb3adjYw3MWM3(3,3),adjYb3Yb3adjYw3Yw3(3,3),adjYb3Yb3adjYw3BMWM3(3,3)             & 
& ,adjYb3Yb3adjYw3TYw3(3,3),adjYb3TYb3adjYb3MBM3(3,3),adjYb3TYb3adjYb3Yb3(3,3)           & 
& ,adjYb3TYb3adjYeYe(3,3),adjYb3TYb3adjYw3MWM3(3,3),adjYb3TYb3adjYw3Yw3(3,3)             & 
& ,adjYdmd2YdadjYd(3,3),adjYdmd2YdadjYu(3,3),adjYdmd2TpYx3CYx3(3,3),adjYdYdmq2adjYd(3,3) & 
& ,adjYdYdmq2adjYu(3,3),adjYdYdadjYdmd2(3,3),adjYdYdadjYdYd(3,3),adjYdYdadjYdTYd(3,3)    

Complex(dp) :: adjYdYdadjYumu2(3,3),adjYdYdadjYuYu(3,3),adjYdYdadjYuTYu(3,3),adjYdTYdadjYdYd(3,3)    & 
& ,adjYdTYdadjYuYu(3,3),adjYdTpYx3CYx3Yd(3,3),adjYdTpYx3CYx3TYd(3,3),adjYdTpTYx3CYx3Yd(3,3)& 
& ,adjYeme2YeadjYb3(3,3),adjYeme2YeadjYe(3,3),adjYeme2YeadjYw3(3,3),adjYeYeml2adjYb3(3,3)& 
& ,adjYeYeml2adjYe(3,3),adjYeYeml2adjYw3(3,3),adjYeYeadjYb3MBM3(3,3),adjYeYeadjYb3mHb32(3,3)& 
& ,adjYeYeadjYb3Yb3(3,3),adjYeYeadjYb3BMBM3(3,3),adjYeYeadjYb3TYb3(3,3),adjYeYeadjYeme2(3,3)& 
& ,adjYeYeadjYeYe(3,3),adjYeYeadjYeTYe(3,3),adjYeYeadjYw3mHw32(3,3),adjYeYeadjYw3MWM3(3,3)& 
& ,adjYeYeadjYw3Yw3(3,3),adjYeYeadjYw3BMWM3(3,3),adjYeYeadjYw3TYw3(3,3),adjYeTYeadjYb3MBM3(3,3)& 
& ,adjYeTYeadjYb3Yb3(3,3),adjYeTYeadjYeYe(3,3),adjYeTYeadjYw3MWM3(3,3),adjYeTYeadjYw3Yw3(3,3)& 
& ,adjYumu2YuadjYd(3,3),adjYumu2YuadjYu(3,3),adjYuYumq2adjYd(3,3),adjYuYumq2adjYu(3,3)   & 
& ,adjYuYuadjYdmd2(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdTYd(3,3),adjYuYuadjYumu2(3,3)    & 
& ,adjYuYuadjYuYu(3,3),adjYuYuadjYuTYu(3,3),adjYuTYuadjYdYd(3,3),adjYuTYuadjYuYu(3,3)    & 
& ,adjYw3mHw32Yw3adjYb3(3,3),adjYw3mHw32Yw3adjYe(3,3),adjYw3mHw32Yw3adjYw3(3,3)          & 
& ,adjYw3Yw3ml2adjYb3(3,3),adjYw3Yw3ml2adjYe(3,3),adjYw3Yw3ml2adjYw3(3,3),               & 
& adjYw3Yw3adjYb3MBM3(3,3),adjYw3Yw3adjYb3mHb32(3,3),adjYw3Yw3adjYb3Yb3(3,3)             & 
& ,adjYw3Yw3adjYb3BMBM3(3,3),adjYw3Yw3adjYb3TYb3(3,3),adjYw3Yw3adjYeme2(3,3)             & 
& ,adjYw3Yw3adjYeYe(3,3),adjYw3Yw3adjYeTYe(3,3),adjYw3Yw3adjYw3mHw32(3,3),               & 
& adjYw3Yw3adjYw3MWM3(3,3),adjYw3Yw3adjYw3Yw3(3,3),adjYw3Yw3adjYw3BMWM3(3,3)             & 
& ,adjYw3Yw3adjYw3TYw3(3,3),adjYw3TYw3adjYb3MBM3(3,3),adjYw3TYw3adjYb3Yb3(3,3)           & 
& ,adjYw3TYw3adjYeYe(3,3),adjYw3TYw3adjYw3MWM3(3,3),adjYw3TYw3adjYw3Yw3(3,3)             & 
& ,adjYx3mHxb32Yx3adjYx3(3,3),adjYx3mHxb32Yx3CYd(3,3),adjYx3Yx3md2adjYx3(3,3)            & 
& ,adjYx3Yx3md2CYd(3,3),adjYx3Yx3adjYx3mHxb32(3,3),adjYx3Yx3adjYx3Yx3(3,3)               & 
& ,adjYx3Yx3adjYx3TYx3(3,3),adjYx3Yx3CYdmq2(3,3),adjYx3TYx3adjYx3Yx3(3,3),               & 
& adjTYb3TYb3TpYb3CYb3(3,3),adjTYb3TYb3TpYw3CYw3(3,3),adjTYdTYdTpYdCYd(3,3)              & 
& ,adjTYeTYeTpYeCYe(3,3),adjTYuTYuTpYuCYu(3,3),adjTYw3TYw3TpYb3CYb3(3,3),adjTYw3TYw3TpYw3CYw3(3,3)& 
& ,adjTYx3TYx3TpYx3CYx3(3,3),CYb3ml2TpYb3CYb3(3,3),CYb3ml2TpYw3CYw3(3,3),CYb3TpYb3mHb32CYb3(3,3)& 
& ,CYb3TpYb3CYb3ml2(3,3),CYb3TpYb3CYb3TpYb3(3,3),CYb3TpYb3CYb3TpTYb3(3,3),               & 
& CYb3TpYeCYeTpYb3(3,3),CYb3TpYeCYeTpTYb3(3,3),CYb3TpYw3mHw32CYw3(3,3),CYb3TpYw3CYw3ml2(3,3)& 
& ,CYb3TpYw3CYw3TpYb3(3,3),CYb3TpYw3CYw3TpTYb3(3,3),CYb3TpTYb3CYb3TpYb3(3,3)             & 
& ,CYb3TpTYeCYeTpYb3(3,3),CYb3TpTYw3CYw3TpYb3(3,3),CYdmq2TpYdadjYx3(3,3),CYdmq2TpYdCYd(3,3)& 
& ,CYdTpYdmd2adjYx3(3,3),CYdTpYdmd2CYd(3,3),CYdTpYdadjYx3mHxb32(3,3),CYdTpYdadjYx3Yx3(3,3)& 
& ,CYdTpYdadjYx3TYx3(3,3),CYdTpYdCYdmq2(3,3),CYdTpYdCYdTpYd(3,3),CYdTpYdCYdTpTYd(3,3)    & 
& ,CYdTpYuCYuTpYd(3,3),CYdTpYuCYuTpTYd(3,3),CYdTpTYdCYdTpYd(3,3),CYdTpTYuCYuTpYd(3,3)    & 
& ,CYeml2TpYeCYe(3,3),CYeTpYeme2CYe(3,3),CYeTpYeCYeml2(3,3),CYumq2TpYuCYu(3,3)           & 
& ,CYuTpYumu2CYu(3,3),CYuTpYuCYumq2(3,3),CYw3ml2TpYb3CYb3(3,3),CYw3ml2TpYw3CYw3(3,3)     & 
& ,CYw3TpYb3mHb32CYb3(3,3),CYw3TpYb3CYb3ml2(3,3),CYw3TpYb3CYb3TpYw3(3,3),CYw3TpYb3CYb3TpTYw3(3,3)& 
& ,CYw3TpYeCYeTpYw3(3,3),CYw3TpYeCYeTpTYw3(3,3),CYw3TpYw3mHw32CYw3(3,3),CYw3TpYw3CYw3ml2(3,3)& 
& ,CYw3TpYw3CYw3TpYw3(3,3),CYw3TpYw3CYw3TpTYw3(3,3),CYw3TpTYb3CYb3TpYw3(3,3)             

Complex(dp) :: CYw3TpTYeCYeTpYw3(3,3),CYw3TpTYw3CYw3TpYw3(3,3),CYx3md2YdadjYd(3,3),CYx3md2TpYx3CYx3(3,3)& 
& ,CYx3YdadjYdTpYx3(3,3),CYx3YdadjYdTpTYx3(3,3),CYx3TYdadjYdTpYx3(3,3),CYx3TpYx3mHxb32CYx3(3,3)& 
& ,CYx3TpYx3CYx3md2(3,3),CYx3TpYx3CYx3Yd(3,3),CYx3TpYx3CYx3TYd(3,3),CYx3TpYx3CYx3TpYx3(3,3)& 
& ,CYx3TpYx3CYx3TpTYx3(3,3),CYx3TpTYx3CYx3Yd(3,3),CYx3TpTYx3CYx3TpYx3(3,3)               & 
& ,CTYdTpYdadjYx3TYx3(3,3),TYb3adjYb3Yb3adjTYb3(3,3),TYb3adjYeYeadjTYb3(3,3)             & 
& ,TYb3adjYw3Yw3adjTYb3(3,3),TYb3TpYb3CTYb3adjYb3(3,3),TYb3TpYeCTYeadjYb3(3,3)           & 
& ,TYb3TpYw3CTYw3adjYb3(3,3),TYdadjYdYdadjTYd(3,3),TYdadjYdadjYx3Yx3(3,3),               & 
& TYdadjYdTpYx3CYx3(3,3),TYdadjYdTpYx3CTYx3(3,3),TYdadjYuYuadjTYd(3,3),TYdTpYdCTYdadjYd(3,3)& 
& ,TYdTpYuCTYuadjYd(3,3),TYeadjYb3Yb3adjTYe(3,3),TYeadjYeYeadjTYe(3,3),TYeadjYw3Yw3adjTYe(3,3)& 
& ,TYeTpYb3CTYb3adjYe(3,3),TYeTpYeCTYeadjYe(3,3),TYeTpYw3CTYw3adjYe(3,3),TYuadjYdYdadjTYu(3,3)& 
& ,TYuadjYuYuadjTYu(3,3),TYuTpYdCTYdadjYu(3,3),TYuTpYuCTYuadjYu(3,3),TYw3adjYb3Yb3adjTYw3(3,3)& 
& ,TYw3adjYeYeadjTYw3(3,3),TYw3adjYw3Yw3adjTYw3(3,3),TYw3TpYb3CTYb3adjYw3(3,3)           & 
& ,TYw3TpYeCTYeadjYw3(3,3),TYw3TpYw3CTYw3adjYw3(3,3),TYx3YdadjTYdadjYx3(3,3)             & 
& ,TYx3adjYx3Yx3adjTYx3(3,3),TYx3CYdTpYdadjTYx3(3,3),TYx3TpYx3CTYx3adjYx3(3,3)           & 
& ,TpYb3CYb3TpYb3CYb3(3,3),TpYb3CYb3TpYw3CYw3(3,3),TpYb3CYb3TpTYb3CTYb3(3,3)             & 
& ,TpYb3CYb3TpTYw3CTYw3(3,3),TpYb3CTYb3TpTYb3CYb3(3,3),TpYb3CTYb3TpTYw3CYw3(3,3)         & 
& ,TpYdadjYdTpTYx3CTYx3(3,3),TpYdadjYx3mHxb32Yx3(3,3),TpYdadjYx3Yx3CYd(3,3)              & 
& ,TpYdadjYx3TYx3CTYd(3,3),TpYdCYdTpYdCYd(3,3),TpYdCYdTpTYdCTYd(3,3),TpYdCTYdTpTYdCYd(3,3)& 
& ,TpYeCYeTpYeCYe(3,3),TpYeCYeTpTYeCTYe(3,3),TpYeCTYeTpTYeCYe(3,3),TpYuCYuTpYuCYu(3,3)   & 
& ,TpYuCYuTpTYuCTYu(3,3),TpYuCTYuTpTYuCYu(3,3),TpYw3CYw3TpYb3CYb3(3,3),TpYw3CYw3TpYw3CYw3(3,3)& 
& ,TpYw3CYw3TpTYb3CTYb3(3,3),TpYw3CYw3TpTYw3CTYw3(3,3),TpYw3CTYw3TpTYb3CYb3(3,3)         & 
& ,TpYw3CTYw3TpTYw3CYw3(3,3),TpYx3CYx3YdadjYd(3,3),TpYx3CYx3TYdadjYd(3,3),               & 
& TpYx3CYx3TpYx3CYx3(3,3),TpYx3CYx3TpTYx3CTYx3(3,3),TpYx3CTYx3TpTYx3CYx3(3,3)            & 
& ,TpTYb3CYb3TpYb3CTYb3(3,3),TpTYb3CYb3TpYw3CTYw3(3,3),TpTYdadjYdTpYx3CTYx3(3,3)         & 
& ,TpTYdadjYx3Yx3CTYd(3,3),TpTYdCYdTpYdCTYd(3,3),TpTYeCYeTpYeCTYe(3,3),TpTYuCYuTpYuCTYu(3,3)& 
& ,TpTYw3CYw3TpYb3CTYb3(3,3),TpTYw3CYw3TpYw3CTYw3(3,3),TpTYx3CYx3TpYx3CTYx3(3,3)         & 
& ,MBM3CYb3TpYb3CYb3TpYb3(3,3),MBM3CYb3TpYb3CYb3TpTYb3(3,3),MBM3CYb3TpYeCYeTpYb3(3,3)    & 
& ,MBM3CYb3TpYeCYeTpTYb3(3,3),MBM3CYb3TpYw3CYw3TpYb3(3,3),MBM3CYb3TpYw3CYw3TpTYb3(3,3)   & 
& ,MBM3CYb3TpTYb3CYb3TpYb3(3,3),MBM3CYb3TpTYeCYeTpYb3(3,3),MBM3CYb3TpTYw3CYw3TpYb3(3,3)  & 
& ,md2YdadjYdYdadjYd(3,3),md2YdadjYdTpYx3CYx3(3,3),md2YdadjYuYuadjYd(3,3),               & 
& md2TpYx3CYx3YdadjYd(3,3),md2TpYx3CYx3TpYx3CYx3(3,3),me2YeadjYb3Yb3adjYe(3,3)           & 
& ,me2YeadjYeYeadjYe(3,3),me2YeadjYw3Yw3adjYe(3,3),mHb32Yb3adjYb3Yb3adjYb3(3,3)          & 
& ,mHb32Yb3adjYeYeadjYb3(3,3),mHb32Yb3adjYw3Yw3adjYb3(3,3),mHw32Yw3adjYb3Yb3adjYw3(3,3)  & 
& ,mHw32Yw3adjYeYeadjYw3(3,3),mHw32Yw3adjYw3Yw3adjYw3(3,3),mHxb32Yx3adjYx3Yx3adjYx3(3,3) & 
& ,mHxb32Yx3CYdTpYdadjYx3(3,3),ml2TpYb3CYb3TpYb3CYb3(3,3),ml2TpYb3CYb3TpYw3CYw3(3,3)     & 
& ,ml2TpYeCYeTpYeCYe(3,3),ml2TpYw3CYw3TpYb3CYb3(3,3),ml2TpYw3CYw3TpYw3CYw3(3,3)          & 
& ,mq2adjYdTpYx3CYx3Yd(3,3),mq2TpYdCYdTpYdCYd(3,3),mq2TpYuCYuTpYuCYu(3,3)               

Complex(dp) :: mu2YuadjYdYdadjYu(3,3),mu2YuadjYuYuadjYu(3,3),MWM3CYw3TpYb3CYb3TpYw3(3,3)              & 
& ,MWM3CYw3TpYb3CYb3TpTYw3(3,3),MWM3CYw3TpYeCYeTpYw3(3,3),MWM3CYw3TpYeCYeTpTYw3(3,3)     & 
& ,MWM3CYw3TpYw3CYw3TpYw3(3,3),MWM3CYw3TpYw3CYw3TpTYw3(3,3),MWM3CYw3TpTYb3CYb3TpYw3(3,3) & 
& ,MWM3CYw3TpTYeCYeTpYw3(3,3),MWM3CYw3TpTYw3CYw3TpYw3(3,3),MXM3CYx3YdadjYdTpYx3(3,3)     & 
& ,MXM3CYx3YdadjYdTpTYx3(3,3),MXM3CYx3TYdadjYdTpYx3(3,3),MXM3CYx3TpYx3CYx3TpYx3(3,3)     & 
& ,MXM3CYx3TpYx3CYx3TpTYx3(3,3),MXM3CYx3TpTYx3CYx3TpYx3(3,3),Yb3ml2adjYb3Yb3adjYb3(3,3)  & 
& ,Yb3ml2adjYeYeadjYb3(3,3),Yb3ml2adjYw3Yw3adjYb3(3,3),Yb3adjYb3mHb32Yb3adjYb3(3,3)      & 
& ,Yb3adjYb3Yb3ml2adjYb3(3,3),Yb3adjYb3Yb3adjYb3MBM3(3,3),Yb3adjYb3Yb3adjYb3mHb32(3,3)   & 
& ,Yb3adjYb3Yb3adjYb3Yb3(3,3),Yb3adjYb3Yb3adjYb3BMBM3(3,3),Yb3adjYb3Yb3adjYb3TYb3(3,3)   & 
& ,Yb3adjYb3Yb3adjYw3Yw3(3,3),Yb3adjYb3Yb3adjYw3TYw3(3,3),Yb3adjYb3TYb3adjYb3MBM3(3,3)   & 
& ,Yb3adjYb3TYb3adjYb3Yb3(3,3),Yb3adjYb3TYb3adjYw3Yw3(3,3),Yb3adjYeme2YeadjYb3(3,3)      & 
& ,Yb3adjYeYeml2adjYb3(3,3),Yb3adjYeYeadjYb3MBM3(3,3),Yb3adjYeYeadjYb3mHb32(3,3)         & 
& ,Yb3adjYeYeadjYb3Yb3(3,3),Yb3adjYeYeadjYb3BMBM3(3,3),Yb3adjYeYeadjYb3TYb3(3,3)         & 
& ,Yb3adjYeYeadjYeYe(3,3),Yb3adjYeYeadjYeTYe(3,3),Yb3adjYeYeadjYw3TYw3(3,3)              & 
& ,Yb3adjYeTYeadjYb3MBM3(3,3),Yb3adjYeTYeadjYb3Yb3(3,3),Yb3adjYeTYeadjYeYe(3,3)          & 
& ,Yb3adjYeTYeadjYw3Yw3(3,3),Yb3adjYw3mHw32Yw3adjYb3(3,3),Yb3adjYw3Yw3ml2adjYb3(3,3)     & 
& ,Yb3adjYw3Yw3adjYb3MBM3(3,3),Yb3adjYw3Yw3adjYb3mHb32(3,3),Yb3adjYw3Yw3adjYb3Yb3(3,3)   & 
& ,Yb3adjYw3Yw3adjYb3BMBM3(3,3),Yb3adjYw3Yw3adjYb3TYb3(3,3),Yb3adjYw3Yw3adjYw3Yw3(3,3)   & 
& ,Yb3adjYw3Yw3adjYw3TYw3(3,3),Yb3adjYw3TYw3adjYb3MBM3(3,3),Yb3adjYw3TYw3adjYb3Yb3(3,3)  & 
& ,Yb3adjYw3TYw3adjYw3Yw3(3,3),Ydmq2adjYdYdadjYd(3,3),Ydmq2adjYdTpYx3CYx3(3,3)           & 
& ,Ydmq2adjYuYuadjYd(3,3),Ydmq2adjYx3Yx3CYd(3,3),YdadjYdmd2YdadjYd(3,3),YdadjYdmd2TpYx3CYx3(3,3)& 
& ,YdadjYdYdmq2adjYd(3,3),YdadjYdYdadjYdmd2(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdTYd(3,3)& 
& ,YdadjYdTYdadjYdYd(3,3),YdadjYdTpYx3CYx3Yd(3,3),YdadjYdTpYx3CYx3TYd(3,3)               & 
& ,YdadjYdTpTYx3CYx3Yd(3,3),YdadjYumu2YuadjYd(3,3),YdadjYuYumq2adjYd(3,3),               & 
& YdadjYuYuadjYdmd2(3,3),YdadjYuYuadjYdYd(3,3),YdadjYuYuadjYdTYd(3,3),YdadjYuYuadjYuYu(3,3)& 
& ,YdadjYuYuadjYuTYu(3,3),YdadjYuTYuadjYdYd(3,3),YdadjYuTYuadjYuYu(3,3),Yeml2adjYb3Yb3adjYe(3,3)& 
& ,Yeml2adjYeYeadjYe(3,3),Yeml2adjYw3Yw3adjYe(3,3),YeadjYb3mHb32Yb3adjYe(3,3)            & 
& ,YeadjYb3Yb3ml2adjYe(3,3),YeadjYb3Yb3adjYb3Yb3(3,3),YeadjYb3Yb3adjYb3TYb3(3,3)         & 
& ,YeadjYb3Yb3adjYeme2(3,3),YeadjYb3Yb3adjYeYe(3,3),YeadjYb3Yb3adjYeTYe(3,3)             & 
& ,YeadjYb3Yb3adjYw3Yw3(3,3),YeadjYb3Yb3adjYw3TYw3(3,3),YeadjYb3TYb3adjYb3Yb3(3,3)       & 
& ,YeadjYb3TYb3adjYeYe(3,3),YeadjYb3TYb3adjYw3Yw3(3,3),YeadjYeme2YeadjYe(3,3)            & 
& ,YeadjYeYeml2adjYe(3,3),YeadjYeYeadjYeme2(3,3),YeadjYeYeadjYeYe(3,3),YeadjYeYeadjYeTYe(3,3)& 
& ,YeadjYeTYeadjYeYe(3,3),YeadjYw3mHw32Yw3adjYe(3,3),YeadjYw3Yw3ml2adjYe(3,3)            & 
& ,YeadjYw3Yw3adjYb3Yb3(3,3),YeadjYw3Yw3adjYb3TYb3(3,3),YeadjYw3Yw3adjYeme2(3,3)         & 
& ,YeadjYw3Yw3adjYeYe(3,3),YeadjYw3Yw3adjYeTYe(3,3),YeadjYw3Yw3adjYw3Yw3(3,3)            & 
& ,YeadjYw3Yw3adjYw3TYw3(3,3),YeadjYw3TYw3adjYb3Yb3(3,3),YeadjYw3TYw3adjYeYe(3,3)        & 
& ,YeadjYw3TYw3adjYw3Yw3(3,3),Yumq2adjYdYdadjYu(3,3),Yumq2adjYuYuadjYu(3,3)              

Complex(dp) :: YuadjYdmd2YdadjYu(3,3),YuadjYdYdmq2adjYu(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYdTYd(3,3)& 
& ,YuadjYdYdadjYumu2(3,3),YuadjYdYdadjYuYu(3,3),YuadjYdYdadjYuTYu(3,3),YuadjYdTYdadjYdYd(3,3)& 
& ,YuadjYdTYdadjYuYu(3,3),YuadjYdTpYx3CYx3Yd(3,3),YuadjYdTpYx3CYx3TYd(3,3)               & 
& ,YuadjYdTpTYx3CYx3Yd(3,3),YuadjYumu2YuadjYu(3,3),YuadjYuYumq2adjYu(3,3),               & 
& YuadjYuYuadjYumu2(3,3),YuadjYuYuadjYuYu(3,3),YuadjYuYuadjYuTYu(3,3),YuadjYuTYuadjYuYu(3,3)& 
& ,Yw3ml2adjYb3Yb3adjYw3(3,3),Yw3ml2adjYeYeadjYw3(3,3),Yw3ml2adjYw3Yw3adjYw3(3,3)        & 
& ,Yw3adjYb3mHb32Yb3adjYw3(3,3),Yw3adjYb3Yb3ml2adjYw3(3,3),Yw3adjYb3Yb3adjYb3Yb3(3,3)    & 
& ,Yw3adjYb3Yb3adjYb3TYb3(3,3),Yw3adjYb3Yb3adjYw3mHw32(3,3),Yw3adjYb3Yb3adjYw3MWM3(3,3)  & 
& ,Yw3adjYb3Yb3adjYw3Yw3(3,3),Yw3adjYb3Yb3adjYw3BMWM3(3,3),Yw3adjYb3Yb3adjYw3TYw3(3,3)   & 
& ,Yw3adjYb3TYb3adjYb3Yb3(3,3),Yw3adjYb3TYb3adjYw3MWM3(3,3),Yw3adjYb3TYb3adjYw3Yw3(3,3)  & 
& ,Yw3adjYeme2YeadjYw3(3,3),Yw3adjYeYeml2adjYw3(3,3),Yw3adjYeYeadjYb3TYb3(3,3)           & 
& ,Yw3adjYeYeadjYeYe(3,3),Yw3adjYeYeadjYeTYe(3,3),Yw3adjYeYeadjYw3mHw32(3,3)             & 
& ,Yw3adjYeYeadjYw3MWM3(3,3),Yw3adjYeYeadjYw3Yw3(3,3),Yw3adjYeYeadjYw3BMWM3(3,3)         & 
& ,Yw3adjYeYeadjYw3TYw3(3,3),Yw3adjYeTYeadjYb3Yb3(3,3),Yw3adjYeTYeadjYeYe(3,3)           & 
& ,Yw3adjYeTYeadjYw3MWM3(3,3),Yw3adjYeTYeadjYw3Yw3(3,3),Yw3adjYw3mHw32Yw3adjYw3(3,3)     & 
& ,Yw3adjYw3Yw3ml2adjYw3(3,3),Yw3adjYw3Yw3adjYb3Yb3(3,3),Yw3adjYw3Yw3adjYb3TYb3(3,3)     & 
& ,Yw3adjYw3Yw3adjYw3mHw32(3,3),Yw3adjYw3Yw3adjYw3MWM3(3,3),Yw3adjYw3Yw3adjYw3Yw3(3,3)   & 
& ,Yw3adjYw3Yw3adjYw3BMWM3(3,3),Yw3adjYw3Yw3adjYw3TYw3(3,3),Yw3adjYw3TYw3adjYb3Yb3(3,3)  & 
& ,Yw3adjYw3TYw3adjYw3MWM3(3,3),Yw3adjYw3TYw3adjYw3Yw3(3,3),Yx3md2adjYx3Yx3adjYx3(3,3)   & 
& ,Yx3md2CYdTpYdadjYx3(3,3),Yx3adjYx3mHxb32Yx3adjYx3(3,3),Yx3adjYx3Yx3md2adjYx3(3,3)     & 
& ,Yx3adjYx3Yx3adjYx3mHxb32(3,3),Yx3adjYx3Yx3adjYx3Yx3(3,3),Yx3adjYx3Yx3adjYx3TYx3(3,3)  & 
& ,Yx3adjYx3TYx3adjYx3Yx3(3,3),Yx3CYdmq2TpYdadjYx3(3,3),Yx3CYdTpYdmd2adjYx3(3,3)         & 
& ,Yx3CYdTpYdadjYx3mHxb32(3,3),Yx3CYdTpYdadjYx3Yx3(3,3),Yx3CYdTpYdadjYx3TYx3(3,3)        & 
& ,Yx3CYdTpYdCYdTpYd(3,3),Yx3CYdTpYdCYdTpTYd(3,3),Yx3CYdTpYuCYuTpYd(3,3),Yx3CYdTpYuCYuTpTYd(3,3)& 
& ,Yx3CYdTpTYdCYdTpYd(3,3),Yx3CYdTpTYuCYuTpYd(3,3),Yx3TYdadjYdadjYx3Yx3(3,3)             & 
& ,BMBM3CYb3TpYb3CYb3TpYb3(3,3),BMBM3CYb3TpYeCYeTpYb3(3,3),BMBM3CYb3TpYw3CYw3TpYb3(3,3)  & 
& ,BMWM3CYw3TpYb3CYb3TpYw3(3,3),BMWM3CYw3TpYeCYeTpYw3(3,3),BMWM3CYw3TpYw3CYw3TpYw3(3,3)  & 
& ,BMXM3CYx3YdadjYdTpYx3(3,3),BMXM3CYx3TpYx3CYx3TpYx3(3,3),CYdTpYdadjYx3mHxb32Yx3(3,3)   & 
& ,TYb3adjYb3Yb3adjYb3MBM3(3,3),TYb3adjYb3Yb3adjYb3Yb3(3,3),TYb3adjYb3Yb3adjYw3Yw3(3,3)  & 
& ,TYb3adjYeYeadjYb3MBM3(3,3),TYb3adjYeYeadjYb3Yb3(3,3),TYb3adjYeYeadjYeYe(3,3)          & 
& ,TYb3adjYeYeadjYw3Yw3(3,3),TYb3adjYw3Yw3adjYb3MBM3(3,3),TYb3adjYw3Yw3adjYb3Yb3(3,3)    & 
& ,TYb3adjYw3Yw3adjYw3Yw3(3,3),TYdadjYdYdadjYdYd(3,3),TYdadjYdTpYx3CYx3Yd(3,3)           & 
& ,TYdadjYuYuadjYdYd(3,3),TYdadjYuYuadjYuYu(3,3),TYeadjYb3Yb3adjYb3Yb3(3,3)              & 
& ,TYeadjYb3Yb3adjYeYe(3,3),TYeadjYb3Yb3adjYw3Yw3(3,3),TYeadjYeYeadjYeYe(3,3)            & 
& ,TYeadjYw3Yw3adjYb3Yb3(3,3),TYeadjYw3Yw3adjYeYe(3,3),TYeadjYw3Yw3adjYw3Yw3(3,3)        & 
& ,TYuadjYdYdadjYdYd(3,3),TYuadjYdYdadjYuYu(3,3),TYuadjYdTpYx3CYx3Yd(3,3),               & 
& TYuadjYuYuadjYuYu(3,3),TYw3adjYb3Yb3adjYb3Yb3(3,3),TYw3adjYb3Yb3adjYw3MWM3(3,3)        

Complex(dp) :: TYw3adjYb3Yb3adjYw3Yw3(3,3),TYw3adjYeYeadjYb3Yb3(3,3),TYw3adjYeYeadjYeYe(3,3)         & 
& ,TYw3adjYeYeadjYw3MWM3(3,3),TYw3adjYeYeadjYw3Yw3(3,3),TYw3adjYw3Yw3adjYb3Yb3(3,3)      & 
& ,TYw3adjYw3Yw3adjYw3MWM3(3,3),TYw3adjYw3Yw3adjYw3Yw3(3,3),TYx3adjYx3Yx3adjYx3Yx3(3,3)  & 
& ,TYx3CYdTpYdadjYx3Yx3(3,3),TYx3CYdTpYdCYdTpYd(3,3),TYx3CYdTpYuCYuTpYd(3,3)             & 
& ,TpYb3mHb32CYb3TpYb3CYb3(3,3),TpYb3mHb32CYb3TpYw3CYw3(3,3),TpYb3CYb3ml2TpYb3CYb3(3,3)  & 
& ,TpYb3CYb3ml2TpYw3CYw3(3,3),TpYb3CYb3TpYb3mHb32CYb3(3,3),TpYb3CYb3TpYb3CYb3ml2(3,3)    & 
& ,TpYb3CYb3TpYw3mHw32CYw3(3,3),TpYb3CYb3TpYw3CYw3ml2(3,3),TpYdmd2adjYx3Yx3CYd(3,3)      & 
& ,TpYdmd2CYdTpYdCYd(3,3),TpYdadjYx3mHxb32Yx3CYd(3,3),TpYdadjYx3Yx3md2CYd(3,3)           & 
& ,TpYdadjYx3Yx3CYdmq2(3,3),TpYdCYdmq2TpYdCYd(3,3),TpYdCYdTpYdmd2CYd(3,3),               & 
& TpYdCYdTpYdCYdmq2(3,3),TpYeme2CYeTpYeCYe(3,3),TpYeCYeml2TpYeCYe(3,3),TpYeCYeTpYeme2CYe(3,3)& 
& ,TpYeCYeTpYeCYeml2(3,3),TpYumu2CYuTpYuCYu(3,3),TpYuCYumq2TpYuCYu(3,3),TpYuCYuTpYumu2CYu(3,3)& 
& ,TpYuCYuTpYuCYumq2(3,3),TpYw3mHw32CYw3TpYb3CYb3(3,3),TpYw3mHw32CYw3TpYw3CYw3(3,3)      & 
& ,TpYw3CYw3ml2TpYb3CYb3(3,3),TpYw3CYw3ml2TpYw3CYw3(3,3),TpYw3CYw3TpYb3mHb32CYb3(3,3)    & 
& ,TpYw3CYw3TpYb3CYb3ml2(3,3),TpYw3CYw3TpYw3mHw32CYw3(3,3),TpYw3CYw3TpYw3CYw3ml2(3,3)    & 
& ,TpYx3mHxb32CYx3YdadjYd(3,3),TpYx3mHxb32CYx3TpYx3CYx3(3,3),TpYx3CYx3md2YdadjYd(3,3)    & 
& ,TpYx3CYx3md2TpYx3CYx3(3,3),TpYx3CYx3TpYx3mHxb32CYx3(3,3),TpYx3CYx3TpYx3CYx3md2(3,3)   & 
& ,TpYx3CYx3TpYx3CYx3Yd(3,3),TpYx3CYx3TpYx3CYx3TYd(3,3),TpYx3CYx3TpTYx3CYx3Yd(3,3)       & 
& ,TpTYx3CYx3TpYx3CYx3Yd(3,3)

Complex(dp) :: Trmd2,Trme2,TrmHx32,TrmHxb32,Trml2,Trmq2,Trmu2,TrYb3adjYb3,TrYdadjYd,TrYeadjYe,       & 
& TrYuadjYu,TrYw3adjYw3,TrYx3adjYx3,TrTYb3adjYb3,TrTYdadjYd,TrTYeadjYe,TrTYuadjYu,       & 
& TrTYw3adjYw3,TrTYx3adjYx3,TrTpTYb3CTYb3,TrTpTYdCTYd,TrTpTYeCTYe,TrTpTYuCTYu,           & 
& TrTpTYw3CTYw3,TrTpTYx3CTYx3,Trmd2YdadjYd,Trme2YeadjYe,TrmHb32Yb3adjYb3,TrmHw32Yw3adjYw3,& 
& TrmHxb32Yx3adjYx3,Trmu2YuadjYu,TrYb3ml2adjYb3,TrYdmq2adjYd,TrYeml2adjYe,               & 
& TrYumq2adjYu,TrYw3ml2adjYw3,TrYx3md2adjYx3

Complex(dp) :: TrmHg32,TrmHw32,TrTpYb3CTYb3,TrTpYdCTYd,TrTpYeCTYe,TrTpYuCTYu,TrTpYw3CTYw3,           & 
& TrTpYx3CTYx3,TrYb3adjYb3Yb3adjYb3,TrYb3adjYb3TYb3adjYb3,TrYb3adjYb3TYb3adjTYb3,        & 
& TrYb3adjYeYeadjYb3,TrYb3adjYeTYeadjYb3,TrYb3adjYeTYeadjTYb3,TrYb3adjYw3Yw3adjYb3,      & 
& TrYb3adjYw3TYw3adjYb3,TrYb3adjYw3TYw3adjTYb3,TrYb3TpTYb3CTYb3adjYb3,TrYb3TpTYeCTYeadjYb3,& 
& TrYb3TpTYw3CTYw3adjYb3,TrYdadjYdYdadjYd,TrYdadjYdTYdadjYd,TrYdadjYdTYdadjTYd,          & 
& TrYdadjYdTpYx3CYx3,TrYdadjYdTpTYx3CTYx3,TrYdadjYuYuadjYd,TrYdadjYuTYuadjYd,            & 
& TrYdadjYuTYuadjTYd,TrYdTpTYdCTYdadjYd,TrYdTpTYuCTYuadjYd,TrYeadjYb3Yb3adjYe,           & 
& TrYeadjYb3TYb3adjYe,TrYeadjYb3TYb3adjTYe,TrYeadjYeYeadjYe,TrYeadjYeTYeadjYe,           & 
& TrYeadjYeTYeadjTYe,TrYeadjYw3Yw3adjYe,TrYeadjYw3TYw3adjYe,TrYeadjYw3TYw3adjTYe,        & 
& TrYeTpTYb3CTYb3adjYe,TrYeTpTYeCTYeadjYe,TrYeTpTYw3CTYw3adjYe,TrYuadjYdYdadjYu,         & 
& TrYuadjYdTYdadjYu,TrYuadjYdTYdadjTYu,TrYuadjYuYuadjYu,TrYuadjYuTYuadjYu,               & 
& TrYuadjYuTYuadjTYu,TrYuTpTYdCTYdadjYu,TrYuTpTYuCTYuadjYu,TrYw3adjYb3Yb3adjYw3,         & 
& TrYw3adjYb3TYb3adjYw3,TrYw3adjYb3TYb3adjTYw3,TrYw3adjYeYeadjYw3,TrYw3adjYeTYeadjYw3,   & 
& TrYw3adjYeTYeadjTYw3,TrYw3adjYw3Yw3adjYw3,TrYw3adjYw3TYw3adjYw3,TrYw3adjYw3TYw3adjTYw3,& 
& TrYw3TpTYb3CTYb3adjYw3,TrYw3TpTYeCTYeadjYw3,TrYw3TpTYw3CTYw3adjYw3,TrYx3adjYx3Yx3adjYx3,& 
& TrYx3adjYx3TYx3adjYx3,TrYx3adjYx3TYx3adjTYx3,TrYx3TpTYx3CTYx3adjYx3,TrCYdTpYdadjYx3TYx3,& 
& TrCTYdTpYdadjYx3TYx3,TrTYb3adjYb3Yb3adjTYb3,TrTYb3adjYeYeadjTYb3,TrTYb3adjYw3Yw3adjTYb3,& 
& TrTYdadjYdYdadjTYd,TrTYdadjYdTpYx3CYx3,TrTYdadjYdTpYx3CTYx3,TrTYdadjYuYuadjTYd,        & 
& TrTYeadjYb3Yb3adjTYe,TrTYeadjYeYeadjTYe,TrTYeadjYw3Yw3adjTYe,TrTYuadjYdYdadjTYu,       & 
& TrTYuadjYuYuadjTYu,TrTYw3adjYb3Yb3adjTYw3,TrTYw3adjYeYeadjTYw3,TrTYw3adjYw3Yw3adjTYw3, & 
& TrTYx3adjYx3Yx3adjTYx3,TrTpYdadjYdTpTYx3CTYx3,TrTpYdadjYx3TYx3CTYd,TrTpYx3CYx3YdadjYd, & 
& TrTpYx3CYx3TYdadjYd,TrTpTYdadjYdTpYx3CTYx3,TrTpTYdadjYx3Yx3CTYd,Trmd2YdadjYdTpYx3CYx3, & 
& Trmd2TpYx3CYx3YdadjYd,Trmq2adjYdTpYx3CYx3Yd,TrYb3ml2adjYb3Yb3adjYb3,TrYb3ml2adjYeYeadjYb3,& 
& TrYb3ml2adjYw3Yw3adjYb3,TrYb3adjYb3mHb32Yb3adjYb3,TrYb3adjYb3Yb3ml2adjYb3,             & 
& TrYb3adjYeme2YeadjYb3,TrYb3adjYeYeml2adjYb3,TrYb3adjYw3mHw32Yw3adjYb3,TrYb3adjYw3Yw3ml2adjYb3,& 
& TrYdmq2adjYdYdadjYd,TrYdmq2adjYdTpYx3CYx3,TrYdmq2adjYuYuadjYd,TrYdadjYdmd2YdadjYd,     & 
& TrYdadjYdmd2TpYx3CYx3,TrYdadjYdYdmq2adjYd,TrYdadjYumu2YuadjYd,TrYdadjYuYumq2adjYd,     & 
& TrYeml2adjYb3Yb3adjYe,TrYeml2adjYeYeadjYe,TrYeml2adjYw3Yw3adjYe,TrYeadjYb3mHb32Yb3adjYe,& 
& TrYeadjYb3Yb3ml2adjYe,TrYeadjYeme2YeadjYe,TrYeadjYeYeml2adjYe,TrYeadjYw3mHw32Yw3adjYe, & 
& TrYeadjYw3Yw3ml2adjYe,TrYumq2adjYdYdadjYu,TrYumq2adjYuYuadjYu,TrYuadjYdmd2YdadjYu,     & 
& TrYuadjYdYdmq2adjYu,TrYuadjYumu2YuadjYu,TrYuadjYuYumq2adjYu,TrYw3ml2adjYb3Yb3adjYw3,   & 
& TrYw3ml2adjYeYeadjYw3,TrYw3ml2adjYw3Yw3adjYw3,TrYw3adjYb3mHb32Yb3adjYw3,               & 
& TrYw3adjYb3Yb3ml2adjYw3,TrYw3adjYeme2YeadjYw3,TrYw3adjYeYeml2adjYw3,TrYw3adjYw3mHw32Yw3adjYw3,& 
& TrYw3adjYw3Yw3ml2adjYw3,TrYx3md2adjYx3Yx3adjYx3,TrYx3adjYx3mHxb32Yx3adjYx3,            & 
& TrYx3adjYx3Yx3md2adjYx3,TrCYdTpYdadjYx3mHxb32Yx3,TrTpYx3mHxb32CYx3YdadjYd,             & 
& TrTpYx3CYx3md2YdadjYd
Complex(dp) :: MnuL(3,3), DMnuL(3,3), betaMnuL1(3,3), betaMNuL2(3,3), gammaL1(3,3), &
 & gammaL2(3,3), gammaHd1, gammaHd2

Real(dp) :: g1p2,g1p3,g2p2,g2p3,g3p2,g3p3

Real(dp) :: g1p4,g1p5,g2p4,g2p5,g3p4,g3p5

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
Call Adjungate(BMXM3,adjBMXM3)
Call Adjungate(BMWM3,adjBMWM3)
Call Adjungate(BMGM3,adjBMGM3)
Call Adjungate(BMBM3,adjBMBM3)
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
 md2adjYx3 = Matmul2(md2,adjYx3,OnlyDiagonal) 
 md2CYd = Matmul2(md2,Conjg(Yd),OnlyDiagonal) 
 me2CYe = Matmul2(me2,Conjg(Ye),OnlyDiagonal) 
 mHb32CYb3 = Matmul2(mHb32,Conjg(Yb3),OnlyDiagonal) 
 mHw32CYw3 = Matmul2(mHw32,Conjg(Yw3),OnlyDiagonal) 
 mHxb32CYx3 = Matmul2(mHxb32,Conjg(Yx3),OnlyDiagonal) 
 ml2adjYb3 = Matmul2(ml2,adjYb3,OnlyDiagonal) 
 ml2adjYe = Matmul2(ml2,adjYe,OnlyDiagonal) 
 ml2adjYw3 = Matmul2(ml2,adjYw3,OnlyDiagonal) 
 mq2adjYd = Matmul2(mq2,adjYd,OnlyDiagonal) 
 mq2adjYu = Matmul2(mq2,adjYu,OnlyDiagonal) 
 mu2CYu = Matmul2(mu2,Conjg(Yu),OnlyDiagonal) 
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
 adjYb3MBM3 = Matmul2(adjYb3,MBM3,OnlyDiagonal) 
 adjYb3mHb32 = Matmul2(adjYb3,mHb32,OnlyDiagonal) 
 adjYb3Yb3 = Matmul2(adjYb3,Yb3,OnlyDiagonal) 
Forall(i2=1:3)  adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3(i2,i2),dp) 
 adjYb3BMBM3 = Matmul2(adjYb3,BMBM3,OnlyDiagonal) 
 adjYb3TYb3 = Matmul2(adjYb3,TYb3,OnlyDiagonal) 
 adjYdmd2 = Matmul2(adjYd,md2,OnlyDiagonal) 
 adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYd(i2,i2) =  Real(adjYdYd(i2,i2),dp) 
 adjYdTYd = Matmul2(adjYd,TYd,OnlyDiagonal) 
 adjYeme2 = Matmul2(adjYe,me2,OnlyDiagonal) 
 adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYe(i2,i2) =  Real(adjYeYe(i2,i2),dp) 
 adjYeTYe = Matmul2(adjYe,TYe,OnlyDiagonal) 
 adjYumu2 = Matmul2(adjYu,mu2,OnlyDiagonal) 
 adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYu(i2,i2) =  Real(adjYuYu(i2,i2),dp) 
 adjYuTYu = Matmul2(adjYu,TYu,OnlyDiagonal) 
 adjYw3mHw32 = Matmul2(adjYw3,mHw32,OnlyDiagonal) 
 adjYw3MWM3 = Matmul2(adjYw3,MWM3,OnlyDiagonal) 
 adjYw3Yw3 = Matmul2(adjYw3,Yw3,OnlyDiagonal) 
Forall(i2=1:3)  adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3(i2,i2),dp) 
 adjYw3BMWM3 = Matmul2(adjYw3,BMWM3,OnlyDiagonal) 
 adjYw3TYw3 = Matmul2(adjYw3,TYw3,OnlyDiagonal) 
 adjYx3mHxb32 = Matmul2(adjYx3,mHxb32,OnlyDiagonal) 
 adjYx3Yx3 = Matmul2(adjYx3,Yx3,OnlyDiagonal) 
Forall(i2=1:3)  adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3(i2,i2),dp) 
 adjYx3TYx3 = Matmul2(adjYx3,TYx3,OnlyDiagonal) 
 CYb3ml2 = Matmul2(Conjg(Yb3),ml2,OnlyDiagonal) 
 CYb3TpYb3 = Matmul2(Conjg(Yb3),Transpose(Yb3),OnlyDiagonal) 
Forall(i2=1:3)  CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3(i2,i2),dp) 
 CYb3TpTYb3 = Matmul2(Conjg(Yb3),Transpose(TYb3),OnlyDiagonal) 
 CYdmq2 = Matmul2(Conjg(Yd),mq2,OnlyDiagonal) 
 CYdTpYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYd(i2,i2) =  Real(CYdTpYd(i2,i2),dp) 
 CYdTpTYd = Matmul2(Conjg(Yd),Transpose(TYd),OnlyDiagonal) 
 CYeml2 = Matmul2(Conjg(Ye),ml2,OnlyDiagonal) 
 CYumq2 = Matmul2(Conjg(Yu),mq2,OnlyDiagonal) 
 CYw3ml2 = Matmul2(Conjg(Yw3),ml2,OnlyDiagonal) 
 CYw3TpYw3 = Matmul2(Conjg(Yw3),Transpose(Yw3),OnlyDiagonal) 
Forall(i2=1:3)  CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3(i2,i2),dp) 
 CYw3TpTYw3 = Matmul2(Conjg(Yw3),Transpose(TYw3),OnlyDiagonal) 
 CYx3md2 = Matmul2(Conjg(Yx3),md2,OnlyDiagonal) 
 CYx3Yd = Matmul2(Conjg(Yx3),Yd,OnlyDiagonal) 
 CYx3TYd = Matmul2(Conjg(Yx3),TYd,OnlyDiagonal) 
 CYx3TpYx3 = Matmul2(Conjg(Yx3),Transpose(Yx3),OnlyDiagonal) 
Forall(i2=1:3)  CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3(i2,i2),dp) 
 CYx3TpTYx3 = Matmul2(Conjg(Yx3),Transpose(TYx3),OnlyDiagonal) 
 TYb3adjYb3 = Matmul2(TYb3,adjYb3,OnlyDiagonal) 
 TYb3adjTYb3 = Matmul2(TYb3,adjTYb3,OnlyDiagonal) 
 TYdadjYd = Matmul2(TYd,adjYd,OnlyDiagonal) 
 TYdadjTYd = Matmul2(TYd,adjTYd,OnlyDiagonal) 
 TYeadjYe = Matmul2(TYe,adjYe,OnlyDiagonal) 
 TYeadjTYe = Matmul2(TYe,adjTYe,OnlyDiagonal) 
 TYuadjYu = Matmul2(TYu,adjYu,OnlyDiagonal) 
 TYuadjTYu = Matmul2(TYu,adjTYu,OnlyDiagonal) 
 TYw3adjYw3 = Matmul2(TYw3,adjYw3,OnlyDiagonal) 
 TYw3adjTYw3 = Matmul2(TYw3,adjTYw3,OnlyDiagonal) 
 TYx3adjYx3 = Matmul2(TYx3,adjYx3,OnlyDiagonal) 
 TYx3adjTYx3 = Matmul2(TYx3,adjTYx3,OnlyDiagonal) 
 TpYb3CYb3 = Matmul2(Transpose(Yb3),Conjg(Yb3),OnlyDiagonal) 
Forall(i2=1:3)  TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3(i2,i2),dp) 
 TpYdCYd = Matmul2(Transpose(Yd),Conjg(Yd),OnlyDiagonal) 
Forall(i2=1:3)  TpYdCYd(i2,i2) =  Real(TpYdCYd(i2,i2),dp) 
 TpYeCYe = Matmul2(Transpose(Ye),Conjg(Ye),OnlyDiagonal) 
Forall(i2=1:3)  TpYeCYe(i2,i2) =  Real(TpYeCYe(i2,i2),dp) 
 TpYuCYu = Matmul2(Transpose(Yu),Conjg(Yu),OnlyDiagonal) 
Forall(i2=1:3)  TpYuCYu(i2,i2) =  Real(TpYuCYu(i2,i2),dp) 
 TpYw3CYw3 = Matmul2(Transpose(Yw3),Conjg(Yw3),OnlyDiagonal) 
Forall(i2=1:3)  TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3(i2,i2),dp) 
 TpYx3CYx3 = Matmul2(Transpose(Yx3),Conjg(Yx3),OnlyDiagonal) 
Forall(i2=1:3)  TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3(i2,i2),dp) 
 TpTYb3CTYb3 = Matmul2(Transpose(TYb3),Conjg(TYb3),OnlyDiagonal) 
 TpTYdCTYd = Matmul2(Transpose(TYd),Conjg(TYd),OnlyDiagonal) 
 TpTYeCTYe = Matmul2(Transpose(TYe),Conjg(TYe),OnlyDiagonal) 
 TpTYuCTYu = Matmul2(Transpose(TYu),Conjg(TYu),OnlyDiagonal) 
 TpTYw3CTYw3 = Matmul2(Transpose(TYw3),Conjg(TYw3),OnlyDiagonal) 
 TpTYx3CTYx3 = Matmul2(Transpose(TYx3),Conjg(TYx3),OnlyDiagonal) 
 MBM3CYb3TpYb3 = Matmul2(MBM3,CYb3TpYb3,OnlyDiagonal) 
 MBM3CYb3TpTYb3 = Matmul2(MBM3,CYb3TpTYb3,OnlyDiagonal) 
 md2YdadjYd = Matmul2(md2,YdadjYd,OnlyDiagonal) 
 md2TpYx3CYx3 = Matmul2(md2,TpYx3CYx3,OnlyDiagonal) 
 me2YeadjYe = Matmul2(me2,YeadjYe,OnlyDiagonal) 
 mHb32Yb3adjYb3 = Matmul2(mHb32,Yb3adjYb3,OnlyDiagonal) 
 mHw32Yw3adjYw3 = Matmul2(mHw32,Yw3adjYw3,OnlyDiagonal) 
 mHxb32Yx3adjYx3 = Matmul2(mHxb32,Yx3adjYx3,OnlyDiagonal) 
 ml2TpYb3CYb3 = Matmul2(ml2,TpYb3CYb3,OnlyDiagonal) 
 ml2TpYeCYe = Matmul2(ml2,TpYeCYe,OnlyDiagonal) 
 ml2TpYw3CYw3 = Matmul2(ml2,TpYw3CYw3,OnlyDiagonal) 
 mq2TpYdCYd = Matmul2(mq2,TpYdCYd,OnlyDiagonal) 
 mq2TpYuCYu = Matmul2(mq2,TpYuCYu,OnlyDiagonal) 
 mu2YuadjYu = Matmul2(mu2,YuadjYu,OnlyDiagonal) 
 MWM3CYw3TpYw3 = Matmul2(MWM3,CYw3TpYw3,OnlyDiagonal) 
 MWM3CYw3TpTYw3 = Matmul2(MWM3,CYw3TpTYw3,OnlyDiagonal) 
 MXM3CYx3TpYx3 = Matmul2(MXM3,CYx3TpYx3,OnlyDiagonal) 
 MXM3CYx3TpTYx3 = Matmul2(MXM3,CYx3TpTYx3,OnlyDiagonal) 
 Yb3ml2adjYb3 = Matmul2(Yb3,ml2adjYb3,OnlyDiagonal) 
 Yb3adjYb3MBM3 = Matmul2(Yb3,adjYb3MBM3,OnlyDiagonal) 
 Yb3adjYb3mHb32 = Matmul2(Yb3,adjYb3mHb32,OnlyDiagonal) 
 Yb3adjYb3Yb3 = Matmul2(Yb3,adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYb3BMBM3 = Matmul2(Yb3,adjYb3BMBM3,OnlyDiagonal) 
 Yb3adjYb3TYb3 = Matmul2(Yb3,adjYb3TYb3,OnlyDiagonal) 
 Yb3adjYeYe = Matmul2(Yb3,adjYeYe,OnlyDiagonal) 
 Yb3adjYeTYe = Matmul2(Yb3,adjYeTYe,OnlyDiagonal) 
 Yb3adjYw3Yw3 = Matmul2(Yb3,adjYw3Yw3,OnlyDiagonal) 
 Yb3adjYw3TYw3 = Matmul2(Yb3,adjYw3TYw3,OnlyDiagonal) 
 Ydmq2adjYd = Matmul2(Yd,mq2adjYd,OnlyDiagonal) 
 YdadjYdmd2 = Matmul2(Yd,adjYdmd2,OnlyDiagonal) 
 YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal) 
 YdadjYdTYd = Matmul2(Yd,adjYdTYd,OnlyDiagonal) 
 YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal) 
 YdadjYuTYu = Matmul2(Yd,adjYuTYu,OnlyDiagonal) 
 Yeml2adjYe = Matmul2(Ye,ml2adjYe,OnlyDiagonal) 
 YeadjYb3Yb3 = Matmul2(Ye,adjYb3Yb3,OnlyDiagonal) 
 YeadjYb3TYb3 = Matmul2(Ye,adjYb3TYb3,OnlyDiagonal) 
 YeadjYeme2 = Matmul2(Ye,adjYeme2,OnlyDiagonal) 
 YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal) 
 YeadjYeTYe = Matmul2(Ye,adjYeTYe,OnlyDiagonal) 
 YeadjYw3Yw3 = Matmul2(Ye,adjYw3Yw3,OnlyDiagonal) 
 YeadjYw3TYw3 = Matmul2(Ye,adjYw3TYw3,OnlyDiagonal) 
 Yumq2adjYu = Matmul2(Yu,mq2adjYu,OnlyDiagonal) 
 YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal) 
 YuadjYdTYd = Matmul2(Yu,adjYdTYd,OnlyDiagonal) 
 YuadjYumu2 = Matmul2(Yu,adjYumu2,OnlyDiagonal) 
 YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal) 
 YuadjYuTYu = Matmul2(Yu,adjYuTYu,OnlyDiagonal) 
 Yw3ml2adjYw3 = Matmul2(Yw3,ml2adjYw3,OnlyDiagonal) 
 Yw3adjYb3Yb3 = Matmul2(Yw3,adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYb3TYb3 = Matmul2(Yw3,adjYb3TYb3,OnlyDiagonal) 
 Yw3adjYeYe = Matmul2(Yw3,adjYeYe,OnlyDiagonal) 
 Yw3adjYeTYe = Matmul2(Yw3,adjYeTYe,OnlyDiagonal) 
 Yw3adjYw3mHw32 = Matmul2(Yw3,adjYw3mHw32,OnlyDiagonal) 
 Yw3adjYw3MWM3 = Matmul2(Yw3,adjYw3MWM3,OnlyDiagonal) 
 Yw3adjYw3Yw3 = Matmul2(Yw3,adjYw3Yw3,OnlyDiagonal) 
 Yw3adjYw3BMWM3 = Matmul2(Yw3,adjYw3BMWM3,OnlyDiagonal) 
 Yw3adjYw3TYw3 = Matmul2(Yw3,adjYw3TYw3,OnlyDiagonal) 
 Yx3md2adjYx3 = Matmul2(Yx3,md2adjYx3,OnlyDiagonal) 
 Yx3adjYx3mHxb32 = Matmul2(Yx3,adjYx3mHxb32,OnlyDiagonal) 
 Yx3adjYx3Yx3 = Matmul2(Yx3,adjYx3Yx3,OnlyDiagonal) 
 Yx3adjYx3TYx3 = Matmul2(Yx3,adjYx3TYx3,OnlyDiagonal) 
 Yx3CYdTpYd = Matmul2(Yx3,CYdTpYd,OnlyDiagonal) 
 Yx3CYdTpTYd = Matmul2(Yx3,CYdTpTYd,OnlyDiagonal) 
 BMBM3CYb3TpYb3 = Matmul2(BMBM3,CYb3TpYb3,OnlyDiagonal) 
 BMWM3CYw3TpYw3 = Matmul2(BMWM3,CYw3TpYw3,OnlyDiagonal) 
 BMXM3CYx3TpYx3 = Matmul2(BMXM3,CYx3TpYx3,OnlyDiagonal) 
 TYb3adjYb3MBM3 = Matmul2(TYb3,adjYb3MBM3,OnlyDiagonal) 
 TYb3adjYb3Yb3 = Matmul2(TYb3,adjYb3Yb3,OnlyDiagonal) 
 TYb3adjYeYe = Matmul2(TYb3,adjYeYe,OnlyDiagonal) 
 TYb3adjYw3Yw3 = Matmul2(TYb3,adjYw3Yw3,OnlyDiagonal) 
 TYdadjYdYd = Matmul2(TYd,adjYdYd,OnlyDiagonal) 
 TYdadjYuYu = Matmul2(TYd,adjYuYu,OnlyDiagonal) 
 TYeadjYb3Yb3 = Matmul2(TYe,adjYb3Yb3,OnlyDiagonal) 
 TYeadjYeYe = Matmul2(TYe,adjYeYe,OnlyDiagonal) 
 TYeadjYw3Yw3 = Matmul2(TYe,adjYw3Yw3,OnlyDiagonal) 
 TYuadjYdYd = Matmul2(TYu,adjYdYd,OnlyDiagonal) 
 TYuadjYuYu = Matmul2(TYu,adjYuYu,OnlyDiagonal) 
 TYw3adjYb3Yb3 = Matmul2(TYw3,adjYb3Yb3,OnlyDiagonal) 
 TYw3adjYeYe = Matmul2(TYw3,adjYeYe,OnlyDiagonal) 
 TYw3adjYw3MWM3 = Matmul2(TYw3,adjYw3MWM3,OnlyDiagonal) 
 TYw3adjYw3Yw3 = Matmul2(TYw3,adjYw3Yw3,OnlyDiagonal) 
 TYx3adjYx3Yx3 = Matmul2(TYx3,adjYx3Yx3,OnlyDiagonal) 
 TYx3CYdTpYd = Matmul2(TYx3,CYdTpYd,OnlyDiagonal) 
 TpYb3mHb32CYb3 = Matmul2(Transpose(Yb3),mHb32CYb3,OnlyDiagonal) 
 TpYb3CYb3ml2 = Matmul2(Transpose(Yb3),CYb3ml2,OnlyDiagonal) 
 TpYdmd2CYd = Matmul2(Transpose(Yd),md2CYd,OnlyDiagonal) 
 TpYdCYdmq2 = Matmul2(Transpose(Yd),CYdmq2,OnlyDiagonal) 
 TpYeme2CYe = Matmul2(Transpose(Ye),me2CYe,OnlyDiagonal) 
 TpYeCYeml2 = Matmul2(Transpose(Ye),CYeml2,OnlyDiagonal) 
 TpYumu2CYu = Matmul2(Transpose(Yu),mu2CYu,OnlyDiagonal) 
 TpYuCYumq2 = Matmul2(Transpose(Yu),CYumq2,OnlyDiagonal) 
 TpYw3mHw32CYw3 = Matmul2(Transpose(Yw3),mHw32CYw3,OnlyDiagonal) 
 TpYw3CYw3ml2 = Matmul2(Transpose(Yw3),CYw3ml2,OnlyDiagonal) 
 TpYx3mHxb32CYx3 = Matmul2(Transpose(Yx3),mHxb32CYx3,OnlyDiagonal) 
 TpYx3CYx3md2 = Matmul2(Transpose(Yx3),CYx3md2,OnlyDiagonal) 
 TpYx3CYx3Yd = Matmul2(Transpose(Yx3),CYx3Yd,OnlyDiagonal) 
 TpYx3CYx3TYd = Matmul2(Transpose(Yx3),CYx3TYd,OnlyDiagonal) 
 TpTYx3CYx3Yd = Matmul2(Transpose(TYx3),CYx3Yd,OnlyDiagonal) 
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
 TrTYb3adjYb3 = Real(cTrace(TYb3adjYb3),dp) 
 TrTYdadjYd = Real(cTrace(TYdadjYd),dp) 
 TrTYeadjYe = Real(cTrace(TYeadjYe),dp) 
 TrTYuadjYu = Real(cTrace(TYuadjYu),dp) 
 TrTYw3adjYw3 = Real(cTrace(TYw3adjYw3),dp) 
 TrTYx3adjYx3 = Real(cTrace(TYx3adjYx3),dp) 
 TrTpTYb3CTYb3 = Real(cTrace(TpTYb3CTYb3),dp) 
 TrTpTYdCTYd = Real(cTrace(TpTYdCTYd),dp) 
 TrTpTYeCTYe = Real(cTrace(TpTYeCTYe),dp) 
 TrTpTYuCTYu = Real(cTrace(TpTYuCTYu),dp) 
 TrTpTYw3CTYw3 = Real(cTrace(TpTYw3CTYw3),dp) 
 TrTpTYx3CTYx3 = Real(cTrace(TpTYx3CTYx3),dp) 
 Trmd2YdadjYd = Real(cTrace(md2YdadjYd),dp) 
 Trme2YeadjYe = Real(cTrace(me2YeadjYe),dp) 
 TrmHb32Yb3adjYb3 = Real(cTrace(mHb32Yb3adjYb3),dp) 
 TrmHw32Yw3adjYw3 = Real(cTrace(mHw32Yw3adjYw3),dp) 
 TrmHxb32Yx3adjYx3 = Real(cTrace(mHxb32Yx3adjYx3),dp) 
 Trmu2YuadjYu = Real(cTrace(mu2YuadjYu),dp) 
 TrYb3ml2adjYb3 = Real(cTrace(Yb3ml2adjYb3),dp) 
 TrYdmq2adjYd = Real(cTrace(Ydmq2adjYd),dp) 
 TrYeml2adjYe = Real(cTrace(Yeml2adjYe),dp) 
 TrYumq2adjYu = Real(cTrace(Yumq2adjYu),dp) 
 TrYw3ml2adjYw3 = Real(cTrace(Yw3ml2adjYw3),dp) 
 TrYx3md2adjYx3 = Real(cTrace(Yx3md2adjYx3),dp) 
 g1p2 =g1**2 
 g1p3 =g1**3 
 g2p2 =g2**2 
 g2p3 =g2**3 
 g3p2 =g3**2 
 g3p3 =g3**3 


If (TwoLoopRGE) Then 
 mHxb32Yx3 = Matmul2(mHxb32,Yx3,OnlyDiagonal) 
 Yb3adjYe = Matmul2(Yb3,adjYe,OnlyDiagonal) 
 Yb3adjYw3 = Matmul2(Yb3,adjYw3,OnlyDiagonal) 
 Yb3adjTYb3 = Matmul2(Yb3,adjTYb3,OnlyDiagonal) 
 Yb3adjTYe = Matmul2(Yb3,adjTYe,OnlyDiagonal) 
 Yb3adjTYw3 = Matmul2(Yb3,adjTYw3,OnlyDiagonal) 
 YdadjYu = Matmul2(Yd,adjYu,OnlyDiagonal) 
 YdadjTYd = Matmul2(Yd,adjTYd,OnlyDiagonal) 
 YdadjTYu = Matmul2(Yd,adjTYu,OnlyDiagonal) 
 YeadjYb3 = Matmul2(Ye,adjYb3,OnlyDiagonal) 
 YeadjYw3 = Matmul2(Ye,adjYw3,OnlyDiagonal) 
 YeadjTYb3 = Matmul2(Ye,adjTYb3,OnlyDiagonal) 
 YeadjTYe = Matmul2(Ye,adjTYe,OnlyDiagonal) 
 YeadjTYw3 = Matmul2(Ye,adjTYw3,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 YuadjTYd = Matmul2(Yu,adjTYd,OnlyDiagonal) 
 YuadjTYu = Matmul2(Yu,adjTYu,OnlyDiagonal) 
 Yw3adjYb3 = Matmul2(Yw3,adjYb3,OnlyDiagonal) 
 Yw3adjYe = Matmul2(Yw3,adjYe,OnlyDiagonal) 
 Yw3adjTYb3 = Matmul2(Yw3,adjTYb3,OnlyDiagonal) 
 Yw3adjTYe = Matmul2(Yw3,adjTYe,OnlyDiagonal) 
 Yw3adjTYw3 = Matmul2(Yw3,adjTYw3,OnlyDiagonal) 
 Yx3adjTYx3 = Matmul2(Yx3,adjTYx3,OnlyDiagonal) 
 Yx3CYd = Matmul2(Yx3,Conjg(Yd),OnlyDiagonal) 
 Yx3CTYd = Matmul2(Yx3,Conjg(TYd),OnlyDiagonal) 
 adjYdadjTYx3 = Matmul2(adjYd,adjTYx3,OnlyDiagonal) 
 adjYdTpYx3 = Matmul2(adjYd,Transpose(Yx3),OnlyDiagonal) 
 adjYdTpTYx3 = Matmul2(adjYd,Transpose(TYx3),OnlyDiagonal) 
 adjTYdadjYx3 = Matmul2(adjTYd,adjYx3,OnlyDiagonal) 
 CYb3TpYw3 = Matmul2(Conjg(Yb3),Transpose(Yw3),OnlyDiagonal) 
 CYb3TpTYw3 = Matmul2(Conjg(Yb3),Transpose(TYw3),OnlyDiagonal) 
 CYeTpYb3 = Matmul2(Conjg(Ye),Transpose(Yb3),OnlyDiagonal) 
 CYeTpYw3 = Matmul2(Conjg(Ye),Transpose(Yw3),OnlyDiagonal) 
 CYeTpTYb3 = Matmul2(Conjg(Ye),Transpose(TYb3),OnlyDiagonal) 
 CYeTpTYw3 = Matmul2(Conjg(Ye),Transpose(TYw3),OnlyDiagonal) 
 CYuTpYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 CYuTpTYd = Matmul2(Conjg(Yu),Transpose(TYd),OnlyDiagonal) 
 CYw3TpYb3 = Matmul2(Conjg(Yw3),Transpose(Yb3),OnlyDiagonal) 
 CYw3TpTYb3 = Matmul2(Conjg(Yw3),Transpose(TYb3),OnlyDiagonal) 
 CTYb3adjYb3 = Matmul2(Conjg(TYb3),adjYb3,OnlyDiagonal) 
 CTYb3adjYe = Matmul2(Conjg(TYb3),adjYe,OnlyDiagonal) 
 CTYb3adjYw3 = Matmul2(Conjg(TYb3),adjYw3,OnlyDiagonal) 
 CTYdadjYd = Matmul2(Conjg(TYd),adjYd,OnlyDiagonal) 
 CTYdadjYu = Matmul2(Conjg(TYd),adjYu,OnlyDiagonal) 
 CTYeadjYb3 = Matmul2(Conjg(TYe),adjYb3,OnlyDiagonal) 
 CTYeadjYe = Matmul2(Conjg(TYe),adjYe,OnlyDiagonal) 
 CTYeadjYw3 = Matmul2(Conjg(TYe),adjYw3,OnlyDiagonal) 
 CTYuadjYd = Matmul2(Conjg(TYu),adjYd,OnlyDiagonal) 
 CTYuadjYu = Matmul2(Conjg(TYu),adjYu,OnlyDiagonal) 
 CTYw3adjYb3 = Matmul2(Conjg(TYw3),adjYb3,OnlyDiagonal) 
 CTYw3adjYe = Matmul2(Conjg(TYw3),adjYe,OnlyDiagonal) 
 CTYw3adjYw3 = Matmul2(Conjg(TYw3),adjYw3,OnlyDiagonal) 
 CTYx3adjYx3 = Matmul2(Conjg(TYx3),adjYx3,OnlyDiagonal) 
 TYb3adjYe = Matmul2(TYb3,adjYe,OnlyDiagonal) 
 TYb3adjYw3 = Matmul2(TYb3,adjYw3,OnlyDiagonal) 
 TYb3adjTYe = Matmul2(TYb3,adjTYe,OnlyDiagonal) 
 TYb3adjTYw3 = Matmul2(TYb3,adjTYw3,OnlyDiagonal) 
 TYdadjYu = Matmul2(TYd,adjYu,OnlyDiagonal) 
 TYdadjTYu = Matmul2(TYd,adjTYu,OnlyDiagonal) 
 TYeadjYb3 = Matmul2(TYe,adjYb3,OnlyDiagonal) 
 TYeadjYw3 = Matmul2(TYe,adjYw3,OnlyDiagonal) 
 TYeadjTYb3 = Matmul2(TYe,adjTYb3,OnlyDiagonal) 
 TYeadjTYw3 = Matmul2(TYe,adjTYw3,OnlyDiagonal) 
 TYuadjYd = Matmul2(TYu,adjYd,OnlyDiagonal) 
 TYuadjTYd = Matmul2(TYu,adjTYd,OnlyDiagonal) 
 TYw3adjYb3 = Matmul2(TYw3,adjYb3,OnlyDiagonal) 
 TYw3adjYe = Matmul2(TYw3,adjYe,OnlyDiagonal) 
 TYw3adjTYb3 = Matmul2(TYw3,adjTYb3,OnlyDiagonal) 
 TYw3adjTYe = Matmul2(TYw3,adjTYe,OnlyDiagonal) 
 TYx3CTYd = Matmul2(TYx3,Conjg(TYd),OnlyDiagonal) 
 TpYb3CTYb3 = Matmul2(Transpose(Yb3),Conjg(TYb3),OnlyDiagonal) 
 TpYdadjYx3 = Matmul2(Transpose(Yd),adjYx3,OnlyDiagonal) 
 TpYdadjTYx3 = Matmul2(Transpose(Yd),adjTYx3,OnlyDiagonal) 
 TpYdCTYd = Matmul2(Transpose(Yd),Conjg(TYd),OnlyDiagonal) 
 TpYeCTYe = Matmul2(Transpose(Ye),Conjg(TYe),OnlyDiagonal) 
 TpYuCTYu = Matmul2(Transpose(Yu),Conjg(TYu),OnlyDiagonal) 
 TpYw3CTYw3 = Matmul2(Transpose(Yw3),Conjg(TYw3),OnlyDiagonal) 
 TpYx3CTYx3 = Matmul2(Transpose(Yx3),Conjg(TYx3),OnlyDiagonal) 
 TpTYb3CYb3 = Matmul2(Transpose(TYb3),Conjg(Yb3),OnlyDiagonal) 
 TpTYdadjYx3 = Matmul2(Transpose(TYd),adjYx3,OnlyDiagonal) 
 TpTYdCYd = Matmul2(Transpose(TYd),Conjg(Yd),OnlyDiagonal) 
 TpTYeCYe = Matmul2(Transpose(TYe),Conjg(Ye),OnlyDiagonal) 
 TpTYuCYu = Matmul2(Transpose(TYu),Conjg(Yu),OnlyDiagonal) 
 TpTYw3CYw3 = Matmul2(Transpose(TYw3),Conjg(Yw3),OnlyDiagonal) 
 TpTYx3CYx3 = Matmul2(Transpose(TYx3),Conjg(Yx3),OnlyDiagonal) 
 md2YdadjYu = Matmul2(md2,YdadjYu,OnlyDiagonal) 
 me2YeadjYb3 = Matmul2(me2,YeadjYb3,OnlyDiagonal) 
 me2YeadjYw3 = Matmul2(me2,YeadjYw3,OnlyDiagonal) 
 mHb32Yb3adjYe = Matmul2(mHb32,Yb3adjYe,OnlyDiagonal) 
 mHb32Yb3adjYw3 = Matmul2(mHb32,Yb3adjYw3,OnlyDiagonal) 
 mHw32Yw3adjYb3 = Matmul2(mHw32,Yw3adjYb3,OnlyDiagonal) 
 mHw32Yw3adjYe = Matmul2(mHw32,Yw3adjYe,OnlyDiagonal) 
 mHxb32Yx3CYd = Matmul2(mHxb32,Yx3CYd,OnlyDiagonal) 
 mq2TpYdadjYx3 = Matmul2(mq2,TpYdadjYx3,OnlyDiagonal) 
 mu2YuadjYd = Matmul2(mu2,YuadjYd,OnlyDiagonal) 
 Yb3ml2adjYe = Matmul2(Yb3,ml2adjYe,OnlyDiagonal) 
 Yb3ml2adjYw3 = Matmul2(Yb3,ml2adjYw3,OnlyDiagonal) 
 Yb3adjYeme2 = Matmul2(Yb3,adjYeme2,OnlyDiagonal) 
 Yb3adjYw3mHw32 = Matmul2(Yb3,adjYw3mHw32,OnlyDiagonal) 
 Yb3adjYw3MWM3 = Matmul2(Yb3,adjYw3MWM3,OnlyDiagonal) 
 Yb3adjYw3BMWM3 = Matmul2(Yb3,adjYw3BMWM3,OnlyDiagonal) 
 Ydmq2adjYu = Matmul2(Yd,mq2adjYu,OnlyDiagonal) 
 YdadjYdTpYx3 = Matmul2(Yd,adjYdTpYx3,OnlyDiagonal) 
 YdadjYdTpTYx3 = Matmul2(Yd,adjYdTpTYx3,OnlyDiagonal) 
 YdadjYumu2 = Matmul2(Yd,adjYumu2,OnlyDiagonal) 
 YdadjTYdadjYx3 = Matmul2(Yd,adjTYdadjYx3,OnlyDiagonal) 
 Yeml2adjYb3 = Matmul2(Ye,ml2adjYb3,OnlyDiagonal) 
 Yeml2adjYw3 = Matmul2(Ye,ml2adjYw3,OnlyDiagonal) 
 YeadjYb3MBM3 = Matmul2(Ye,adjYb3MBM3,OnlyDiagonal) 
 YeadjYb3mHb32 = Matmul2(Ye,adjYb3mHb32,OnlyDiagonal) 
 YeadjYb3BMBM3 = Matmul2(Ye,adjYb3BMBM3,OnlyDiagonal) 
 YeadjYw3mHw32 = Matmul2(Ye,adjYw3mHw32,OnlyDiagonal) 
 YeadjYw3MWM3 = Matmul2(Ye,adjYw3MWM3,OnlyDiagonal) 
 YeadjYw3BMWM3 = Matmul2(Ye,adjYw3BMWM3,OnlyDiagonal) 
 Yumq2adjYd = Matmul2(Yu,mq2adjYd,OnlyDiagonal) 
 YuadjYdmd2 = Matmul2(Yu,adjYdmd2,OnlyDiagonal) 
 Yw3ml2adjYb3 = Matmul2(Yw3,ml2adjYb3,OnlyDiagonal) 
 Yw3ml2adjYe = Matmul2(Yw3,ml2adjYe,OnlyDiagonal) 
 Yw3adjYb3MBM3 = Matmul2(Yw3,adjYb3MBM3,OnlyDiagonal) 
 Yw3adjYb3mHb32 = Matmul2(Yw3,adjYb3mHb32,OnlyDiagonal) 
 Yw3adjYb3BMBM3 = Matmul2(Yw3,adjYb3BMBM3,OnlyDiagonal) 
 Yw3adjYeme2 = Matmul2(Yw3,adjYeme2,OnlyDiagonal) 
 Yx3md2CYd = Matmul2(Yx3,md2CYd,OnlyDiagonal) 
 Yx3CYdmq2 = Matmul2(Yx3,CYdmq2,OnlyDiagonal) 
 adjYb3Yb3adjYb3 = Matmul2(adjYb3,Yb3adjYb3,OnlyDiagonal) 
 adjYb3Yb3adjYe = Matmul2(adjYb3,Yb3adjYe,OnlyDiagonal) 
 adjYb3Yb3adjYw3 = Matmul2(adjYb3,Yb3adjYw3,OnlyDiagonal) 
 adjYb3Yb3adjTYb3 = Matmul2(adjYb3,Yb3adjTYb3,OnlyDiagonal) 
 adjYb3Yb3adjTYe = Matmul2(adjYb3,Yb3adjTYe,OnlyDiagonal) 
 adjYb3Yb3adjTYw3 = Matmul2(adjYb3,Yb3adjTYw3,OnlyDiagonal) 
 adjYb3TYb3adjYb3 = Matmul2(adjYb3,TYb3adjYb3,OnlyDiagonal) 
 adjYb3TYb3adjYe = Matmul2(adjYb3,TYb3adjYe,OnlyDiagonal) 
 adjYb3TYb3adjYw3 = Matmul2(adjYb3,TYb3adjYw3,OnlyDiagonal) 
 adjYb3TYb3adjTYb3 = Matmul2(adjYb3,TYb3adjTYb3,OnlyDiagonal) 
 adjYb3TYb3adjTYe = Matmul2(adjYb3,TYb3adjTYe,OnlyDiagonal) 
 adjYb3TYb3adjTYw3 = Matmul2(adjYb3,TYb3adjTYw3,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdYdadjYu = Matmul2(adjYd,YdadjYu,OnlyDiagonal) 
 adjYdYdadjTYd = Matmul2(adjYd,YdadjTYd,OnlyDiagonal) 
 adjYdYdadjTYu = Matmul2(adjYd,YdadjTYu,OnlyDiagonal) 
 adjYdadjYx3Yx3 = Matmul2(adjYd,adjYx3Yx3,OnlyDiagonal) 
 adjYdTYdadjYd = Matmul2(adjYd,TYdadjYd,OnlyDiagonal) 
 adjYdTYdadjYu = Matmul2(adjYd,TYdadjYu,OnlyDiagonal) 
 adjYdTYdadjTYd = Matmul2(adjYd,TYdadjTYd,OnlyDiagonal) 
 adjYdTYdadjTYu = Matmul2(adjYd,TYdadjTYu,OnlyDiagonal) 
 adjYdTpYx3CYx3 = Matmul2(adjYd,TpYx3CYx3,OnlyDiagonal) 
 adjYdTpYx3CTYx3 = Matmul2(adjYd,TpYx3CTYx3,OnlyDiagonal) 
 adjYdTpTYx3CTYx3 = Matmul2(adjYd,TpTYx3CTYx3,OnlyDiagonal) 
 adjYeYeadjYb3 = Matmul2(adjYe,YeadjYb3,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYeYeadjYw3 = Matmul2(adjYe,YeadjYw3,OnlyDiagonal) 
 adjYeYeadjTYb3 = Matmul2(adjYe,YeadjTYb3,OnlyDiagonal) 
 adjYeYeadjTYe = Matmul2(adjYe,YeadjTYe,OnlyDiagonal) 
 adjYeYeadjTYw3 = Matmul2(adjYe,YeadjTYw3,OnlyDiagonal) 
 adjYeTYeadjYb3 = Matmul2(adjYe,TYeadjYb3,OnlyDiagonal) 
 adjYeTYeadjYe = Matmul2(adjYe,TYeadjYe,OnlyDiagonal) 
 adjYeTYeadjYw3 = Matmul2(adjYe,TYeadjYw3,OnlyDiagonal) 
 adjYeTYeadjTYb3 = Matmul2(adjYe,TYeadjTYb3,OnlyDiagonal) 
 adjYeTYeadjTYe = Matmul2(adjYe,TYeadjTYe,OnlyDiagonal) 
 adjYeTYeadjTYw3 = Matmul2(adjYe,TYeadjTYw3,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYuYuadjTYd = Matmul2(adjYu,YuadjTYd,OnlyDiagonal) 
 adjYuYuadjTYu = Matmul2(adjYu,YuadjTYu,OnlyDiagonal) 
 adjYuTYuadjYd = Matmul2(adjYu,TYuadjYd,OnlyDiagonal) 
 adjYuTYuadjYu = Matmul2(adjYu,TYuadjYu,OnlyDiagonal) 
 adjYuTYuadjTYd = Matmul2(adjYu,TYuadjTYd,OnlyDiagonal) 
 adjYuTYuadjTYu = Matmul2(adjYu,TYuadjTYu,OnlyDiagonal) 
 adjYw3Yw3adjYb3 = Matmul2(adjYw3,Yw3adjYb3,OnlyDiagonal) 
 adjYw3Yw3adjYe = Matmul2(adjYw3,Yw3adjYe,OnlyDiagonal) 
 adjYw3Yw3adjYw3 = Matmul2(adjYw3,Yw3adjYw3,OnlyDiagonal) 
 adjYw3Yw3adjTYb3 = Matmul2(adjYw3,Yw3adjTYb3,OnlyDiagonal) 
 adjYw3Yw3adjTYe = Matmul2(adjYw3,Yw3adjTYe,OnlyDiagonal) 
 adjYw3Yw3adjTYw3 = Matmul2(adjYw3,Yw3adjTYw3,OnlyDiagonal) 
 adjYw3TYw3adjYb3 = Matmul2(adjYw3,TYw3adjYb3,OnlyDiagonal) 
 adjYw3TYw3adjYe = Matmul2(adjYw3,TYw3adjYe,OnlyDiagonal) 
 adjYw3TYw3adjYw3 = Matmul2(adjYw3,TYw3adjYw3,OnlyDiagonal) 
 adjYw3TYw3adjTYb3 = Matmul2(adjYw3,TYw3adjTYb3,OnlyDiagonal) 
 adjYw3TYw3adjTYe = Matmul2(adjYw3,TYw3adjTYe,OnlyDiagonal) 
 adjYw3TYw3adjTYw3 = Matmul2(adjYw3,TYw3adjTYw3,OnlyDiagonal) 
 adjYx3mHxb32Yx3 = Matmul2(adjYx3,mHxb32Yx3,OnlyDiagonal) 
 adjYx3Yx3adjYx3 = Matmul2(adjYx3,Yx3adjYx3,OnlyDiagonal) 
 adjYx3Yx3adjTYx3 = Matmul2(adjYx3,Yx3adjTYx3,OnlyDiagonal) 
 adjYx3Yx3CYd = Matmul2(adjYx3,Yx3CYd,OnlyDiagonal) 
 adjYx3Yx3CTYd = Matmul2(adjYx3,Yx3CTYd,OnlyDiagonal) 
 adjYx3TYx3adjYx3 = Matmul2(adjYx3,TYx3adjYx3,OnlyDiagonal) 
 adjYx3TYx3adjTYx3 = Matmul2(adjYx3,TYx3adjTYx3,OnlyDiagonal) 
 adjYx3TYx3CTYd = Matmul2(adjYx3,TYx3CTYd,OnlyDiagonal) 
 adjTYb3TYb3adjYb3 = Matmul2(adjTYb3,TYb3adjYb3,OnlyDiagonal) 
 adjTYb3TYb3adjYe = Matmul2(adjTYb3,TYb3adjYe,OnlyDiagonal) 
 adjTYb3TYb3adjYw3 = Matmul2(adjTYb3,TYb3adjYw3,OnlyDiagonal) 
 adjTYdTYdadjYd = Matmul2(adjTYd,TYdadjYd,OnlyDiagonal) 
 adjTYdTYdadjYu = Matmul2(adjTYd,TYdadjYu,OnlyDiagonal) 
 adjTYeTYeadjYb3 = Matmul2(adjTYe,TYeadjYb3,OnlyDiagonal) 
 adjTYeTYeadjYe = Matmul2(adjTYe,TYeadjYe,OnlyDiagonal) 
 adjTYeTYeadjYw3 = Matmul2(adjTYe,TYeadjYw3,OnlyDiagonal) 
 adjTYuTYuadjYd = Matmul2(adjTYu,TYuadjYd,OnlyDiagonal) 
 adjTYuTYuadjYu = Matmul2(adjTYu,TYuadjYu,OnlyDiagonal) 
 adjTYw3TYw3adjYb3 = Matmul2(adjTYw3,TYw3adjYb3,OnlyDiagonal) 
 adjTYw3TYw3adjYe = Matmul2(adjTYw3,TYw3adjYe,OnlyDiagonal) 
 adjTYw3TYw3adjYw3 = Matmul2(adjTYw3,TYw3adjYw3,OnlyDiagonal) 
 adjTYx3TYx3adjYx3 = Matmul2(adjTYx3,TYx3adjYx3,OnlyDiagonal) 
 CYb3TpYb3CYb3 = Matmul2(Conjg(Yb3),TpYb3CYb3,OnlyDiagonal) 
 CYb3TpYb3CTYb3 = Matmul2(Conjg(Yb3),TpYb3CTYb3,OnlyDiagonal) 
 CYb3TpYw3CYw3 = Matmul2(Conjg(Yb3),TpYw3CYw3,OnlyDiagonal) 
 CYb3TpYw3CTYw3 = Matmul2(Conjg(Yb3),TpYw3CTYw3,OnlyDiagonal) 
 CYb3TpTYb3CTYb3 = Matmul2(Conjg(Yb3),TpTYb3CTYb3,OnlyDiagonal) 
 CYb3TpTYw3CTYw3 = Matmul2(Conjg(Yb3),TpTYw3CTYw3,OnlyDiagonal) 
 CYdTpYdadjYx3 = Matmul2(Conjg(Yd),TpYdadjYx3,OnlyDiagonal) 
 CYdTpYdadjTYx3 = Matmul2(Conjg(Yd),TpYdadjTYx3,OnlyDiagonal) 
 CYdTpYdCYd = Matmul2(Conjg(Yd),TpYdCYd,OnlyDiagonal) 
 CYdTpYdCTYd = Matmul2(Conjg(Yd),TpYdCTYd,OnlyDiagonal) 
 CYdTpTYdCTYd = Matmul2(Conjg(Yd),TpTYdCTYd,OnlyDiagonal) 
 CYeTpYeCYe = Matmul2(Conjg(Ye),TpYeCYe,OnlyDiagonal) 
 CYeTpYeCTYe = Matmul2(Conjg(Ye),TpYeCTYe,OnlyDiagonal) 
 CYeTpTYeCTYe = Matmul2(Conjg(Ye),TpTYeCTYe,OnlyDiagonal) 
 CYuTpYuCYu = Matmul2(Conjg(Yu),TpYuCYu,OnlyDiagonal) 
 CYuTpYuCTYu = Matmul2(Conjg(Yu),TpYuCTYu,OnlyDiagonal) 
 CYuTpTYuCTYu = Matmul2(Conjg(Yu),TpTYuCTYu,OnlyDiagonal) 
 CYw3TpYb3CYb3 = Matmul2(Conjg(Yw3),TpYb3CYb3,OnlyDiagonal) 
 CYw3TpYb3CTYb3 = Matmul2(Conjg(Yw3),TpYb3CTYb3,OnlyDiagonal) 
 CYw3TpYw3CYw3 = Matmul2(Conjg(Yw3),TpYw3CYw3,OnlyDiagonal) 
 CYw3TpYw3CTYw3 = Matmul2(Conjg(Yw3),TpYw3CTYw3,OnlyDiagonal) 
 CYw3TpTYb3CTYb3 = Matmul2(Conjg(Yw3),TpTYb3CTYb3,OnlyDiagonal) 
 CYw3TpTYw3CTYw3 = Matmul2(Conjg(Yw3),TpTYw3CTYw3,OnlyDiagonal) 
 CYx3YdadjYd = Matmul2(Conjg(Yx3),YdadjYd,OnlyDiagonal) 
 CYx3TYdadjYd = Matmul2(Conjg(Yx3),TYdadjYd,OnlyDiagonal) 
 CYx3TpYx3CYx3 = Matmul2(Conjg(Yx3),TpYx3CYx3,OnlyDiagonal) 
 CYx3TpYx3CTYx3 = Matmul2(Conjg(Yx3),TpYx3CTYx3,OnlyDiagonal) 
 CYx3TpTYx3CTYx3 = Matmul2(Conjg(Yx3),TpTYx3CTYx3,OnlyDiagonal) 
 CTYb3TpTYb3CYb3 = Matmul2(Conjg(TYb3),TpTYb3CYb3,OnlyDiagonal) 
 CTYb3TpTYw3CYw3 = Matmul2(Conjg(TYb3),TpTYw3CYw3,OnlyDiagonal) 
 CTYdTpTYdadjYx3 = Matmul2(Conjg(TYd),TpTYdadjYx3,OnlyDiagonal) 
 CTYdTpTYdCYd = Matmul2(Conjg(TYd),TpTYdCYd,OnlyDiagonal) 
 CTYeTpTYeCYe = Matmul2(Conjg(TYe),TpTYeCYe,OnlyDiagonal) 
 CTYuTpTYuCYu = Matmul2(Conjg(TYu),TpTYuCYu,OnlyDiagonal) 
 CTYw3TpTYb3CYb3 = Matmul2(Conjg(TYw3),TpTYb3CYb3,OnlyDiagonal) 
 CTYw3TpTYw3CYw3 = Matmul2(Conjg(TYw3),TpTYw3CYw3,OnlyDiagonal) 
 CTYx3TpTYx3CYx3 = Matmul2(Conjg(TYx3),TpTYx3CYx3,OnlyDiagonal) 
 TYb3adjYw3MWM3 = Matmul2(TYb3,adjYw3MWM3,OnlyDiagonal) 
 TYb3TpYb3CYb3 = Matmul2(TYb3,TpYb3CYb3,OnlyDiagonal) 
 TYb3TpYw3CYw3 = Matmul2(TYb3,TpYw3CYw3,OnlyDiagonal) 
 TYdadjYdadjTYx3 = Matmul2(TYd,adjYdadjTYx3,OnlyDiagonal) 
 TYdadjYdTpYx3 = Matmul2(TYd,adjYdTpYx3,OnlyDiagonal) 
 TYdTpYdCYd = Matmul2(TYd,TpYdCYd,OnlyDiagonal) 
 TYeadjYb3MBM3 = Matmul2(TYe,adjYb3MBM3,OnlyDiagonal) 
 TYeadjYw3MWM3 = Matmul2(TYe,adjYw3MWM3,OnlyDiagonal) 
 TYeTpYeCYe = Matmul2(TYe,TpYeCYe,OnlyDiagonal) 
 TYuTpYuCYu = Matmul2(TYu,TpYuCYu,OnlyDiagonal) 
 TYw3adjYb3MBM3 = Matmul2(TYw3,adjYb3MBM3,OnlyDiagonal) 
 TYw3TpYb3CYb3 = Matmul2(TYw3,TpYb3CYb3,OnlyDiagonal) 
 TYw3TpYw3CYw3 = Matmul2(TYw3,TpYw3CYw3,OnlyDiagonal) 
 TYx3TpYx3CYx3 = Matmul2(TYx3,TpYx3CYx3,OnlyDiagonal) 
 TpYb3CYb3TpYb3 = Matmul2(Transpose(Yb3),CYb3TpYb3,OnlyDiagonal) 
 TpYb3CYb3TpYw3 = Matmul2(Transpose(Yb3),CYb3TpYw3,OnlyDiagonal) 
 TpYb3CYb3TpTYb3 = Matmul2(Transpose(Yb3),CYb3TpTYb3,OnlyDiagonal) 
 TpYb3CYb3TpTYw3 = Matmul2(Transpose(Yb3),CYb3TpTYw3,OnlyDiagonal) 
 TpYb3CTYb3adjYb3 = Matmul2(Transpose(Yb3),CTYb3adjYb3,OnlyDiagonal) 
 TpYb3CTYb3adjYe = Matmul2(Transpose(Yb3),CTYb3adjYe,OnlyDiagonal) 
 TpYb3CTYb3adjYw3 = Matmul2(Transpose(Yb3),CTYb3adjYw3,OnlyDiagonal) 
 TpYdmd2adjYx3 = Matmul2(Transpose(Yd),md2adjYx3,OnlyDiagonal) 
 TpYdadjYx3mHxb32 = Matmul2(Transpose(Yd),adjYx3mHxb32,OnlyDiagonal) 
 TpYdadjYx3Yx3 = Matmul2(Transpose(Yd),adjYx3Yx3,OnlyDiagonal) 
 TpYdadjYx3TYx3 = Matmul2(Transpose(Yd),adjYx3TYx3,OnlyDiagonal) 
 TpYdCYdTpYd = Matmul2(Transpose(Yd),CYdTpYd,OnlyDiagonal) 
 TpYdCYdTpTYd = Matmul2(Transpose(Yd),CYdTpTYd,OnlyDiagonal) 
 TpYdCTYdadjYd = Matmul2(Transpose(Yd),CTYdadjYd,OnlyDiagonal) 
 TpYdCTYdadjYu = Matmul2(Transpose(Yd),CTYdadjYu,OnlyDiagonal) 
 TpYeCYeTpYb3 = Matmul2(Transpose(Ye),CYeTpYb3,OnlyDiagonal) 
 TpYeCYeTpYw3 = Matmul2(Transpose(Ye),CYeTpYw3,OnlyDiagonal) 
 TpYeCYeTpTYb3 = Matmul2(Transpose(Ye),CYeTpTYb3,OnlyDiagonal) 
 TpYeCYeTpTYw3 = Matmul2(Transpose(Ye),CYeTpTYw3,OnlyDiagonal) 
 TpYeCTYeadjYb3 = Matmul2(Transpose(Ye),CTYeadjYb3,OnlyDiagonal) 
 TpYeCTYeadjYe = Matmul2(Transpose(Ye),CTYeadjYe,OnlyDiagonal) 
 TpYeCTYeadjYw3 = Matmul2(Transpose(Ye),CTYeadjYw3,OnlyDiagonal) 
 TpYuCYuTpYd = Matmul2(Transpose(Yu),CYuTpYd,OnlyDiagonal) 
 TpYuCYuTpTYd = Matmul2(Transpose(Yu),CYuTpTYd,OnlyDiagonal) 
 TpYuCTYuadjYd = Matmul2(Transpose(Yu),CTYuadjYd,OnlyDiagonal) 
 TpYuCTYuadjYu = Matmul2(Transpose(Yu),CTYuadjYu,OnlyDiagonal) 
 TpYw3CYw3TpYb3 = Matmul2(Transpose(Yw3),CYw3TpYb3,OnlyDiagonal) 
 TpYw3CYw3TpYw3 = Matmul2(Transpose(Yw3),CYw3TpYw3,OnlyDiagonal) 
 TpYw3CYw3TpTYb3 = Matmul2(Transpose(Yw3),CYw3TpTYb3,OnlyDiagonal) 
 TpYw3CYw3TpTYw3 = Matmul2(Transpose(Yw3),CYw3TpTYw3,OnlyDiagonal) 
 TpYw3CTYw3adjYb3 = Matmul2(Transpose(Yw3),CTYw3adjYb3,OnlyDiagonal) 
 TpYw3CTYw3adjYe = Matmul2(Transpose(Yw3),CTYw3adjYe,OnlyDiagonal) 
 TpYw3CTYw3adjYw3 = Matmul2(Transpose(Yw3),CTYw3adjYw3,OnlyDiagonal) 
 TpYx3CYx3TpYx3 = Matmul2(Transpose(Yx3),CYx3TpYx3,OnlyDiagonal) 
 TpYx3CYx3TpTYx3 = Matmul2(Transpose(Yx3),CYx3TpTYx3,OnlyDiagonal) 
 TpYx3CTYx3adjYx3 = Matmul2(Transpose(Yx3),CTYx3adjYx3,OnlyDiagonal) 
 TpTYb3CYb3TpYb3 = Matmul2(Transpose(TYb3),CYb3TpYb3,OnlyDiagonal) 
 TpTYb3CYb3TpYw3 = Matmul2(Transpose(TYb3),CYb3TpYw3,OnlyDiagonal) 
 TpTYb3CTYb3adjYb3 = Matmul2(Transpose(TYb3),CTYb3adjYb3,OnlyDiagonal) 
 TpTYb3CTYb3adjYe = Matmul2(Transpose(TYb3),CTYb3adjYe,OnlyDiagonal) 
 TpTYb3CTYb3adjYw3 = Matmul2(Transpose(TYb3),CTYb3adjYw3,OnlyDiagonal) 
 TpTYdCYdTpYd = Matmul2(Transpose(TYd),CYdTpYd,OnlyDiagonal) 
 TpTYdCTYdadjYd = Matmul2(Transpose(TYd),CTYdadjYd,OnlyDiagonal) 
 TpTYdCTYdadjYu = Matmul2(Transpose(TYd),CTYdadjYu,OnlyDiagonal) 
 TpTYeCYeTpYb3 = Matmul2(Transpose(TYe),CYeTpYb3,OnlyDiagonal) 
 TpTYeCYeTpYw3 = Matmul2(Transpose(TYe),CYeTpYw3,OnlyDiagonal) 
 TpTYeCTYeadjYb3 = Matmul2(Transpose(TYe),CTYeadjYb3,OnlyDiagonal) 
 TpTYeCTYeadjYe = Matmul2(Transpose(TYe),CTYeadjYe,OnlyDiagonal) 
 TpTYeCTYeadjYw3 = Matmul2(Transpose(TYe),CTYeadjYw3,OnlyDiagonal) 
 TpTYuCYuTpYd = Matmul2(Transpose(TYu),CYuTpYd,OnlyDiagonal) 
 TpTYuCTYuadjYd = Matmul2(Transpose(TYu),CTYuadjYd,OnlyDiagonal) 
 TpTYuCTYuadjYu = Matmul2(Transpose(TYu),CTYuadjYu,OnlyDiagonal) 
 TpTYw3CYw3TpYb3 = Matmul2(Transpose(TYw3),CYw3TpYb3,OnlyDiagonal) 
 TpTYw3CYw3TpYw3 = Matmul2(Transpose(TYw3),CYw3TpYw3,OnlyDiagonal) 
 TpTYw3CTYw3adjYb3 = Matmul2(Transpose(TYw3),CTYw3adjYb3,OnlyDiagonal) 
 TpTYw3CTYw3adjYe = Matmul2(Transpose(TYw3),CTYw3adjYe,OnlyDiagonal) 
 TpTYw3CTYw3adjYw3 = Matmul2(Transpose(TYw3),CTYw3adjYw3,OnlyDiagonal) 
 TpTYx3CYx3TpYx3 = Matmul2(Transpose(TYx3),CYx3TpYx3,OnlyDiagonal) 
 TpTYx3CTYx3adjYx3 = Matmul2(Transpose(TYx3),CTYx3adjYx3,OnlyDiagonal) 
 md2adjYx3Yx3adjYx3 = Matmul2(md2,adjYx3Yx3adjYx3,OnlyDiagonal) 
 md2adjYx3Yx3CYd = Matmul2(md2,adjYx3Yx3CYd,OnlyDiagonal) 
 md2CYdTpYdadjYx3 = Matmul2(md2,CYdTpYdadjYx3,OnlyDiagonal) 
 md2CYdTpYdCYd = Matmul2(md2,CYdTpYdCYd,OnlyDiagonal) 
 me2CYeTpYeCYe = Matmul2(me2,CYeTpYeCYe,OnlyDiagonal) 
 mHb32CYb3TpYb3CYb3 = Matmul2(mHb32,CYb3TpYb3CYb3,OnlyDiagonal) 
 mHb32CYb3TpYw3CYw3 = Matmul2(mHb32,CYb3TpYw3CYw3,OnlyDiagonal) 
 mHw32CYw3TpYb3CYb3 = Matmul2(mHw32,CYw3TpYb3CYb3,OnlyDiagonal) 
 mHw32CYw3TpYw3CYw3 = Matmul2(mHw32,CYw3TpYw3CYw3,OnlyDiagonal) 
 mHxb32CYx3YdadjYd = Matmul2(mHxb32,CYx3YdadjYd,OnlyDiagonal) 
 mHxb32CYx3TpYx3CYx3 = Matmul2(mHxb32,CYx3TpYx3CYx3,OnlyDiagonal) 
 ml2adjYb3Yb3adjYb3 = Matmul2(ml2,adjYb3Yb3adjYb3,OnlyDiagonal) 
 ml2adjYb3Yb3adjYe = Matmul2(ml2,adjYb3Yb3adjYe,OnlyDiagonal) 
 ml2adjYb3Yb3adjYw3 = Matmul2(ml2,adjYb3Yb3adjYw3,OnlyDiagonal) 
 ml2adjYeYeadjYb3 = Matmul2(ml2,adjYeYeadjYb3,OnlyDiagonal) 
 ml2adjYeYeadjYe = Matmul2(ml2,adjYeYeadjYe,OnlyDiagonal) 
 ml2adjYeYeadjYw3 = Matmul2(ml2,adjYeYeadjYw3,OnlyDiagonal) 
 ml2adjYw3Yw3adjYb3 = Matmul2(ml2,adjYw3Yw3adjYb3,OnlyDiagonal) 
 ml2adjYw3Yw3adjYe = Matmul2(ml2,adjYw3Yw3adjYe,OnlyDiagonal) 
 ml2adjYw3Yw3adjYw3 = Matmul2(ml2,adjYw3Yw3adjYw3,OnlyDiagonal) 
 mq2adjYdYdadjYd = Matmul2(mq2,adjYdYdadjYd,OnlyDiagonal) 
 mq2adjYdYdadjYu = Matmul2(mq2,adjYdYdadjYu,OnlyDiagonal) 
 mq2adjYdTpYx3CYx3 = Matmul2(mq2,adjYdTpYx3CYx3,OnlyDiagonal) 
 mq2adjYuYuadjYd = Matmul2(mq2,adjYuYuadjYd,OnlyDiagonal) 
 mq2adjYuYuadjYu = Matmul2(mq2,adjYuYuadjYu,OnlyDiagonal) 
 mq2adjYx3Yx3CYd = Matmul2(mq2,adjYx3Yx3CYd,OnlyDiagonal) 
 mu2CYuTpYuCYu = Matmul2(mu2,CYuTpYuCYu,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3 = Matmul2(Yb3,adjYb3Yb3adjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYb3Yb3adjYb3(i2,i2) =  Real(Yb3adjYb3Yb3adjYb3(i2,i2),dp) 
 Yb3adjYb3TYb3adjYb3 = Matmul2(Yb3,adjYb3TYb3adjYb3,OnlyDiagonal) 
 Yb3adjYb3TYb3adjTYb3 = Matmul2(Yb3,adjYb3TYb3adjTYb3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3 = Matmul2(Yb3,adjYeYeadjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYeYeadjYb3(i2,i2) =  Real(Yb3adjYeYeadjYb3(i2,i2),dp) 
 Yb3adjYeTYeadjYb3 = Matmul2(Yb3,adjYeTYeadjYb3,OnlyDiagonal) 
 Yb3adjYeTYeadjTYb3 = Matmul2(Yb3,adjYeTYeadjTYb3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3 = Matmul2(Yb3,adjYw3Yw3adjYb3,OnlyDiagonal) 
Forall(i2=1:3)  Yb3adjYw3Yw3adjYb3(i2,i2) =  Real(Yb3adjYw3Yw3adjYb3(i2,i2),dp) 
 Yb3adjYw3TYw3adjYb3 = Matmul2(Yb3,adjYw3TYw3adjYb3,OnlyDiagonal) 
 Yb3adjYw3TYw3adjTYb3 = Matmul2(Yb3,adjYw3TYw3adjTYb3,OnlyDiagonal) 
 Yb3adjTYb3TYb3adjYb3 = Matmul2(Yb3,adjTYb3TYb3adjYb3,OnlyDiagonal) 
 Yb3adjTYeTYeadjYb3 = Matmul2(Yb3,adjTYeTYeadjYb3,OnlyDiagonal) 
 Yb3adjTYw3TYw3adjYb3 = Matmul2(Yb3,adjTYw3TYw3adjYb3,OnlyDiagonal) 
 Yb3TpTYb3CTYb3adjYb3 = Matmul2(Yb3,TpTYb3CTYb3adjYb3,OnlyDiagonal) 
 Yb3TpTYeCTYeadjYb3 = Matmul2(Yb3,TpTYeCTYeadjYb3,OnlyDiagonal) 
 Yb3TpTYw3CTYw3adjYb3 = Matmul2(Yb3,TpTYw3CTYw3adjYb3,OnlyDiagonal) 
 YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYdYdadjYd(i2,i2) =  Real(YdadjYdYdadjYd(i2,i2),dp) 
 YdadjYdTYdadjYd = Matmul2(Yd,adjYdTYdadjYd,OnlyDiagonal) 
 YdadjYdTYdadjTYd = Matmul2(Yd,adjYdTYdadjTYd,OnlyDiagonal) 
 YdadjYdTpYx3CYx3 = Matmul2(Yd,adjYdTpYx3CYx3,OnlyDiagonal) 
 YdadjYdTpTYx3CTYx3 = Matmul2(Yd,adjYdTpTYx3CTYx3,OnlyDiagonal) 
 YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal) 
Forall(i2=1:3)  YdadjYuYuadjYd(i2,i2) =  Real(YdadjYuYuadjYd(i2,i2),dp) 
 YdadjYuTYuadjYd = Matmul2(Yd,adjYuTYuadjYd,OnlyDiagonal) 
 YdadjYuTYuadjTYd = Matmul2(Yd,adjYuTYuadjTYd,OnlyDiagonal) 
 YdadjTYdTYdadjYd = Matmul2(Yd,adjTYdTYdadjYd,OnlyDiagonal) 
 YdadjTYuTYuadjYd = Matmul2(Yd,adjTYuTYuadjYd,OnlyDiagonal) 
 YdTpTYdCTYdadjYd = Matmul2(Yd,TpTYdCTYdadjYd,OnlyDiagonal) 
 YdTpTYuCTYuadjYd = Matmul2(Yd,TpTYuCTYuadjYd,OnlyDiagonal) 
 YeadjYb3Yb3adjYe = Matmul2(Ye,adjYb3Yb3adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYb3Yb3adjYe(i2,i2) =  Real(YeadjYb3Yb3adjYe(i2,i2),dp) 
 YeadjYb3TYb3adjYe = Matmul2(Ye,adjYb3TYb3adjYe,OnlyDiagonal) 
 YeadjYb3TYb3adjTYe = Matmul2(Ye,adjYb3TYb3adjTYe,OnlyDiagonal) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYeYeadjYe(i2,i2) =  Real(YeadjYeYeadjYe(i2,i2),dp) 
 YeadjYeTYeadjYe = Matmul2(Ye,adjYeTYeadjYe,OnlyDiagonal) 
 YeadjYeTYeadjTYe = Matmul2(Ye,adjYeTYeadjTYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYe = Matmul2(Ye,adjYw3Yw3adjYe,OnlyDiagonal) 
Forall(i2=1:3)  YeadjYw3Yw3adjYe(i2,i2) =  Real(YeadjYw3Yw3adjYe(i2,i2),dp) 
 YeadjYw3TYw3adjYe = Matmul2(Ye,adjYw3TYw3adjYe,OnlyDiagonal) 
 YeadjYw3TYw3adjTYe = Matmul2(Ye,adjYw3TYw3adjTYe,OnlyDiagonal) 
 YeadjTYb3TYb3adjYe = Matmul2(Ye,adjTYb3TYb3adjYe,OnlyDiagonal) 
 YeadjTYeTYeadjYe = Matmul2(Ye,adjTYeTYeadjYe,OnlyDiagonal) 
 YeadjTYw3TYw3adjYe = Matmul2(Ye,adjTYw3TYw3adjYe,OnlyDiagonal) 
 YeTpTYb3CTYb3adjYe = Matmul2(Ye,TpTYb3CTYb3adjYe,OnlyDiagonal) 
 YeTpTYeCTYeadjYe = Matmul2(Ye,TpTYeCTYeadjYe,OnlyDiagonal) 
 YeTpTYw3CTYw3adjYe = Matmul2(Ye,TpTYw3CTYw3adjYe,OnlyDiagonal) 
 YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYdYdadjYu(i2,i2) =  Real(YuadjYdYdadjYu(i2,i2),dp) 
 YuadjYdTYdadjYu = Matmul2(Yu,adjYdTYdadjYu,OnlyDiagonal) 
 YuadjYdTYdadjTYu = Matmul2(Yu,adjYdTYdadjTYu,OnlyDiagonal) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
Forall(i2=1:3)  YuadjYuYuadjYu(i2,i2) =  Real(YuadjYuYuadjYu(i2,i2),dp) 
 YuadjYuTYuadjYu = Matmul2(Yu,adjYuTYuadjYu,OnlyDiagonal) 
 YuadjYuTYuadjTYu = Matmul2(Yu,adjYuTYuadjTYu,OnlyDiagonal) 
 YuadjTYdTYdadjYu = Matmul2(Yu,adjTYdTYdadjYu,OnlyDiagonal) 
 YuadjTYuTYuadjYu = Matmul2(Yu,adjTYuTYuadjYu,OnlyDiagonal) 
 YuTpTYdCTYdadjYu = Matmul2(Yu,TpTYdCTYdadjYu,OnlyDiagonal) 
 YuTpTYuCTYuadjYu = Matmul2(Yu,TpTYuCTYuadjYu,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3 = Matmul2(Yw3,adjYb3Yb3adjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYb3Yb3adjYw3(i2,i2) =  Real(Yw3adjYb3Yb3adjYw3(i2,i2),dp) 
 Yw3adjYb3TYb3adjYw3 = Matmul2(Yw3,adjYb3TYb3adjYw3,OnlyDiagonal) 
 Yw3adjYb3TYb3adjTYw3 = Matmul2(Yw3,adjYb3TYb3adjTYw3,OnlyDiagonal) 
 Yw3adjYeYeadjYw3 = Matmul2(Yw3,adjYeYeadjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYeYeadjYw3(i2,i2) =  Real(Yw3adjYeYeadjYw3(i2,i2),dp) 
 Yw3adjYeTYeadjYw3 = Matmul2(Yw3,adjYeTYeadjYw3,OnlyDiagonal) 
 Yw3adjYeTYeadjTYw3 = Matmul2(Yw3,adjYeTYeadjTYw3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3 = Matmul2(Yw3,adjYw3Yw3adjYw3,OnlyDiagonal) 
Forall(i2=1:3)  Yw3adjYw3Yw3adjYw3(i2,i2) =  Real(Yw3adjYw3Yw3adjYw3(i2,i2),dp) 
 Yw3adjYw3TYw3adjYw3 = Matmul2(Yw3,adjYw3TYw3adjYw3,OnlyDiagonal) 
 Yw3adjYw3TYw3adjTYw3 = Matmul2(Yw3,adjYw3TYw3adjTYw3,OnlyDiagonal) 
 Yw3adjTYb3TYb3adjYw3 = Matmul2(Yw3,adjTYb3TYb3adjYw3,OnlyDiagonal) 
 Yw3adjTYeTYeadjYw3 = Matmul2(Yw3,adjTYeTYeadjYw3,OnlyDiagonal) 
 Yw3adjTYw3TYw3adjYw3 = Matmul2(Yw3,adjTYw3TYw3adjYw3,OnlyDiagonal) 
 Yw3TpTYb3CTYb3adjYw3 = Matmul2(Yw3,TpTYb3CTYb3adjYw3,OnlyDiagonal) 
 Yw3TpTYeCTYeadjYw3 = Matmul2(Yw3,TpTYeCTYeadjYw3,OnlyDiagonal) 
 Yw3TpTYw3CTYw3adjYw3 = Matmul2(Yw3,TpTYw3CTYw3adjYw3,OnlyDiagonal) 
 Yx3adjYx3Yx3adjYx3 = Matmul2(Yx3,adjYx3Yx3adjYx3,OnlyDiagonal) 
Forall(i2=1:3)  Yx3adjYx3Yx3adjYx3(i2,i2) =  Real(Yx3adjYx3Yx3adjYx3(i2,i2),dp) 
 Yx3adjYx3TYx3adjYx3 = Matmul2(Yx3,adjYx3TYx3adjYx3,OnlyDiagonal) 
 Yx3adjYx3TYx3adjTYx3 = Matmul2(Yx3,adjYx3TYx3adjTYx3,OnlyDiagonal) 
 Yx3adjTYx3TYx3adjYx3 = Matmul2(Yx3,adjTYx3TYx3adjYx3,OnlyDiagonal) 
 Yx3CYdTpYdadjYx3 = Matmul2(Yx3,CYdTpYdadjYx3,OnlyDiagonal) 
Forall(i2=1:3)  Yx3CYdTpYdadjYx3(i2,i2) =  Real(Yx3CYdTpYdadjYx3(i2,i2),dp) 
 Yx3CTYdTpTYdadjYx3 = Matmul2(Yx3,CTYdTpTYdadjYx3,OnlyDiagonal) 
 Yx3TYdadjYdadjTYx3 = Matmul2(Yx3,TYdadjYdadjTYx3,OnlyDiagonal) 
 Yx3TpTYx3CTYx3adjYx3 = Matmul2(Yx3,TpTYx3CTYx3adjYx3,OnlyDiagonal) 
 adjYb3mHb32Yb3adjYb3 = Matmul2(adjYb3,mHb32Yb3adjYb3,OnlyDiagonal) 
 adjYb3mHb32Yb3adjYe = Matmul2(adjYb3,mHb32Yb3adjYe,OnlyDiagonal) 
 adjYb3mHb32Yb3adjYw3 = Matmul2(adjYb3,mHb32Yb3adjYw3,OnlyDiagonal) 
 adjYb3Yb3ml2adjYb3 = Matmul2(adjYb3,Yb3ml2adjYb3,OnlyDiagonal) 
 adjYb3Yb3ml2adjYe = Matmul2(adjYb3,Yb3ml2adjYe,OnlyDiagonal) 
 adjYb3Yb3ml2adjYw3 = Matmul2(adjYb3,Yb3ml2adjYw3,OnlyDiagonal) 
 adjYb3Yb3adjYb3MBM3 = Matmul2(adjYb3,Yb3adjYb3MBM3,OnlyDiagonal) 
 adjYb3Yb3adjYb3mHb32 = Matmul2(adjYb3,Yb3adjYb3mHb32,OnlyDiagonal) 
 adjYb3Yb3adjYb3Yb3 = Matmul2(adjYb3,Yb3adjYb3Yb3,OnlyDiagonal) 
Forall(i2=1:3)  adjYb3Yb3adjYb3Yb3(i2,i2) =  Real(adjYb3Yb3adjYb3Yb3(i2,i2),dp) 
 adjYb3Yb3adjYb3BMBM3 = Matmul2(adjYb3,Yb3adjYb3BMBM3,OnlyDiagonal) 
 adjYb3Yb3adjYb3TYb3 = Matmul2(adjYb3,Yb3adjYb3TYb3,OnlyDiagonal) 
 adjYb3Yb3adjYeme2 = Matmul2(adjYb3,Yb3adjYeme2,OnlyDiagonal) 
 adjYb3Yb3adjYeYe = Matmul2(adjYb3,Yb3adjYeYe,OnlyDiagonal) 
 adjYb3Yb3adjYeTYe = Matmul2(adjYb3,Yb3adjYeTYe,OnlyDiagonal) 
 adjYb3Yb3adjYw3mHw32 = Matmul2(adjYb3,Yb3adjYw3mHw32,OnlyDiagonal) 
 adjYb3Yb3adjYw3MWM3 = Matmul2(adjYb3,Yb3adjYw3MWM3,OnlyDiagonal) 
 adjYb3Yb3adjYw3Yw3 = Matmul2(adjYb3,Yb3adjYw3Yw3,OnlyDiagonal) 
 adjYb3Yb3adjYw3BMWM3 = Matmul2(adjYb3,Yb3adjYw3BMWM3,OnlyDiagonal) 
 adjYb3Yb3adjYw3TYw3 = Matmul2(adjYb3,Yb3adjYw3TYw3,OnlyDiagonal) 
 adjYb3TYb3adjYb3MBM3 = Matmul2(adjYb3,TYb3adjYb3MBM3,OnlyDiagonal) 
 adjYb3TYb3adjYb3Yb3 = Matmul2(adjYb3,TYb3adjYb3Yb3,OnlyDiagonal) 
 adjYb3TYb3adjYeYe = Matmul2(adjYb3,TYb3adjYeYe,OnlyDiagonal) 
 adjYb3TYb3adjYw3MWM3 = Matmul2(adjYb3,TYb3adjYw3MWM3,OnlyDiagonal) 
 adjYb3TYb3adjYw3Yw3 = Matmul2(adjYb3,TYb3adjYw3Yw3,OnlyDiagonal) 
 adjYdmd2YdadjYd = Matmul2(adjYd,md2YdadjYd,OnlyDiagonal) 
 adjYdmd2YdadjYu = Matmul2(adjYd,md2YdadjYu,OnlyDiagonal) 
 adjYdmd2TpYx3CYx3 = Matmul2(adjYd,md2TpYx3CYx3,OnlyDiagonal) 
 adjYdYdmq2adjYd = Matmul2(adjYd,Ydmq2adjYd,OnlyDiagonal) 
 adjYdYdmq2adjYu = Matmul2(adjYd,Ydmq2adjYu,OnlyDiagonal) 
 adjYdYdadjYdmd2 = Matmul2(adjYd,YdadjYdmd2,OnlyDiagonal) 
 adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdYdadjYdYd(i2,i2) =  Real(adjYdYdadjYdYd(i2,i2),dp) 
 adjYdYdadjYdTYd = Matmul2(adjYd,YdadjYdTYd,OnlyDiagonal) 
 adjYdYdadjYumu2 = Matmul2(adjYd,YdadjYumu2,OnlyDiagonal) 
 adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal) 
 adjYdYdadjYuTYu = Matmul2(adjYd,YdadjYuTYu,OnlyDiagonal) 
 adjYdTYdadjYdYd = Matmul2(adjYd,TYdadjYdYd,OnlyDiagonal) 
 adjYdTYdadjYuYu = Matmul2(adjYd,TYdadjYuYu,OnlyDiagonal) 
 adjYdTpYx3CYx3Yd = Matmul2(adjYd,TpYx3CYx3Yd,OnlyDiagonal) 
Forall(i2=1:3)  adjYdTpYx3CYx3Yd(i2,i2) =  Real(adjYdTpYx3CYx3Yd(i2,i2),dp) 
 adjYdTpYx3CYx3TYd = Matmul2(adjYd,TpYx3CYx3TYd,OnlyDiagonal) 
 adjYdTpTYx3CYx3Yd = Matmul2(adjYd,TpTYx3CYx3Yd,OnlyDiagonal) 
 adjYeme2YeadjYb3 = Matmul2(adjYe,me2YeadjYb3,OnlyDiagonal) 
 adjYeme2YeadjYe = Matmul2(adjYe,me2YeadjYe,OnlyDiagonal) 
 adjYeme2YeadjYw3 = Matmul2(adjYe,me2YeadjYw3,OnlyDiagonal) 
 adjYeYeml2adjYb3 = Matmul2(adjYe,Yeml2adjYb3,OnlyDiagonal) 
 adjYeYeml2adjYe = Matmul2(adjYe,Yeml2adjYe,OnlyDiagonal) 
 adjYeYeml2adjYw3 = Matmul2(adjYe,Yeml2adjYw3,OnlyDiagonal) 
 adjYeYeadjYb3MBM3 = Matmul2(adjYe,YeadjYb3MBM3,OnlyDiagonal) 
 adjYeYeadjYb3mHb32 = Matmul2(adjYe,YeadjYb3mHb32,OnlyDiagonal) 
 adjYeYeadjYb3Yb3 = Matmul2(adjYe,YeadjYb3Yb3,OnlyDiagonal) 
 adjYeYeadjYb3BMBM3 = Matmul2(adjYe,YeadjYb3BMBM3,OnlyDiagonal) 
 adjYeYeadjYb3TYb3 = Matmul2(adjYe,YeadjYb3TYb3,OnlyDiagonal) 
 adjYeYeadjYeme2 = Matmul2(adjYe,YeadjYeme2,OnlyDiagonal) 
 adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal) 
Forall(i2=1:3)  adjYeYeadjYeYe(i2,i2) =  Real(adjYeYeadjYeYe(i2,i2),dp) 
 adjYeYeadjYeTYe = Matmul2(adjYe,YeadjYeTYe,OnlyDiagonal) 
 adjYeYeadjYw3mHw32 = Matmul2(adjYe,YeadjYw3mHw32,OnlyDiagonal) 
 adjYeYeadjYw3MWM3 = Matmul2(adjYe,YeadjYw3MWM3,OnlyDiagonal) 
 adjYeYeadjYw3Yw3 = Matmul2(adjYe,YeadjYw3Yw3,OnlyDiagonal) 
 adjYeYeadjYw3BMWM3 = Matmul2(adjYe,YeadjYw3BMWM3,OnlyDiagonal) 
 adjYeYeadjYw3TYw3 = Matmul2(adjYe,YeadjYw3TYw3,OnlyDiagonal) 
 adjYeTYeadjYb3MBM3 = Matmul2(adjYe,TYeadjYb3MBM3,OnlyDiagonal) 
 adjYeTYeadjYb3Yb3 = Matmul2(adjYe,TYeadjYb3Yb3,OnlyDiagonal) 
 adjYeTYeadjYeYe = Matmul2(adjYe,TYeadjYeYe,OnlyDiagonal) 
 adjYeTYeadjYw3MWM3 = Matmul2(adjYe,TYeadjYw3MWM3,OnlyDiagonal) 
 adjYeTYeadjYw3Yw3 = Matmul2(adjYe,TYeadjYw3Yw3,OnlyDiagonal) 
 adjYumu2YuadjYd = Matmul2(adjYu,mu2YuadjYd,OnlyDiagonal) 
 adjYumu2YuadjYu = Matmul2(adjYu,mu2YuadjYu,OnlyDiagonal) 
 adjYuYumq2adjYd = Matmul2(adjYu,Yumq2adjYd,OnlyDiagonal) 
 adjYuYumq2adjYu = Matmul2(adjYu,Yumq2adjYu,OnlyDiagonal) 
 adjYuYuadjYdmd2 = Matmul2(adjYu,YuadjYdmd2,OnlyDiagonal) 
 adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal) 
 adjYuYuadjYdTYd = Matmul2(adjYu,YuadjYdTYd,OnlyDiagonal) 
 adjYuYuadjYumu2 = Matmul2(adjYu,YuadjYumu2,OnlyDiagonal) 
 adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal) 
Forall(i2=1:3)  adjYuYuadjYuYu(i2,i2) =  Real(adjYuYuadjYuYu(i2,i2),dp) 
 adjYuYuadjYuTYu = Matmul2(adjYu,YuadjYuTYu,OnlyDiagonal) 
 adjYuTYuadjYdYd = Matmul2(adjYu,TYuadjYdYd,OnlyDiagonal) 
 adjYuTYuadjYuYu = Matmul2(adjYu,TYuadjYuYu,OnlyDiagonal) 
 adjYw3mHw32Yw3adjYb3 = Matmul2(adjYw3,mHw32Yw3adjYb3,OnlyDiagonal) 
 adjYw3mHw32Yw3adjYe = Matmul2(adjYw3,mHw32Yw3adjYe,OnlyDiagonal) 
 adjYw3mHw32Yw3adjYw3 = Matmul2(adjYw3,mHw32Yw3adjYw3,OnlyDiagonal) 
 adjYw3Yw3ml2adjYb3 = Matmul2(adjYw3,Yw3ml2adjYb3,OnlyDiagonal) 
 adjYw3Yw3ml2adjYe = Matmul2(adjYw3,Yw3ml2adjYe,OnlyDiagonal) 
 adjYw3Yw3ml2adjYw3 = Matmul2(adjYw3,Yw3ml2adjYw3,OnlyDiagonal) 
 adjYw3Yw3adjYb3MBM3 = Matmul2(adjYw3,Yw3adjYb3MBM3,OnlyDiagonal) 
 adjYw3Yw3adjYb3mHb32 = Matmul2(adjYw3,Yw3adjYb3mHb32,OnlyDiagonal) 
 adjYw3Yw3adjYb3Yb3 = Matmul2(adjYw3,Yw3adjYb3Yb3,OnlyDiagonal) 
 adjYw3Yw3adjYb3BMBM3 = Matmul2(adjYw3,Yw3adjYb3BMBM3,OnlyDiagonal) 
 adjYw3Yw3adjYb3TYb3 = Matmul2(adjYw3,Yw3adjYb3TYb3,OnlyDiagonal) 
 adjYw3Yw3adjYeme2 = Matmul2(adjYw3,Yw3adjYeme2,OnlyDiagonal) 
 adjYw3Yw3adjYeYe = Matmul2(adjYw3,Yw3adjYeYe,OnlyDiagonal) 
 adjYw3Yw3adjYeTYe = Matmul2(adjYw3,Yw3adjYeTYe,OnlyDiagonal) 
 adjYw3Yw3adjYw3mHw32 = Matmul2(adjYw3,Yw3adjYw3mHw32,OnlyDiagonal) 
 adjYw3Yw3adjYw3MWM3 = Matmul2(adjYw3,Yw3adjYw3MWM3,OnlyDiagonal) 
 adjYw3Yw3adjYw3Yw3 = Matmul2(adjYw3,Yw3adjYw3Yw3,OnlyDiagonal) 
Forall(i2=1:3)  adjYw3Yw3adjYw3Yw3(i2,i2) =  Real(adjYw3Yw3adjYw3Yw3(i2,i2),dp) 
 adjYw3Yw3adjYw3BMWM3 = Matmul2(adjYw3,Yw3adjYw3BMWM3,OnlyDiagonal) 
 adjYw3Yw3adjYw3TYw3 = Matmul2(adjYw3,Yw3adjYw3TYw3,OnlyDiagonal) 
 adjYw3TYw3adjYb3MBM3 = Matmul2(adjYw3,TYw3adjYb3MBM3,OnlyDiagonal) 
 adjYw3TYw3adjYb3Yb3 = Matmul2(adjYw3,TYw3adjYb3Yb3,OnlyDiagonal) 
 adjYw3TYw3adjYeYe = Matmul2(adjYw3,TYw3adjYeYe,OnlyDiagonal) 
 adjYw3TYw3adjYw3MWM3 = Matmul2(adjYw3,TYw3adjYw3MWM3,OnlyDiagonal) 
 adjYw3TYw3adjYw3Yw3 = Matmul2(adjYw3,TYw3adjYw3Yw3,OnlyDiagonal) 
 adjYx3mHxb32Yx3adjYx3 = Matmul2(adjYx3,mHxb32Yx3adjYx3,OnlyDiagonal) 
 adjYx3mHxb32Yx3CYd = Matmul2(adjYx3,mHxb32Yx3CYd,OnlyDiagonal) 
 adjYx3Yx3md2adjYx3 = Matmul2(adjYx3,Yx3md2adjYx3,OnlyDiagonal) 
 adjYx3Yx3md2CYd = Matmul2(adjYx3,Yx3md2CYd,OnlyDiagonal) 
 adjYx3Yx3adjYx3mHxb32 = Matmul2(adjYx3,Yx3adjYx3mHxb32,OnlyDiagonal) 
 adjYx3Yx3adjYx3Yx3 = Matmul2(adjYx3,Yx3adjYx3Yx3,OnlyDiagonal) 
Forall(i2=1:3)  adjYx3Yx3adjYx3Yx3(i2,i2) =  Real(adjYx3Yx3adjYx3Yx3(i2,i2),dp) 
 adjYx3Yx3adjYx3TYx3 = Matmul2(adjYx3,Yx3adjYx3TYx3,OnlyDiagonal) 
 adjYx3Yx3CYdmq2 = Matmul2(adjYx3,Yx3CYdmq2,OnlyDiagonal) 
 adjYx3TYx3adjYx3Yx3 = Matmul2(adjYx3,TYx3adjYx3Yx3,OnlyDiagonal) 
 adjTYb3TYb3TpYb3CYb3 = Matmul2(adjTYb3,TYb3TpYb3CYb3,OnlyDiagonal) 
 adjTYb3TYb3TpYw3CYw3 = Matmul2(adjTYb3,TYb3TpYw3CYw3,OnlyDiagonal) 
 adjTYdTYdTpYdCYd = Matmul2(adjTYd,TYdTpYdCYd,OnlyDiagonal) 
 adjTYeTYeTpYeCYe = Matmul2(adjTYe,TYeTpYeCYe,OnlyDiagonal) 
 adjTYuTYuTpYuCYu = Matmul2(adjTYu,TYuTpYuCYu,OnlyDiagonal) 
 adjTYw3TYw3TpYb3CYb3 = Matmul2(adjTYw3,TYw3TpYb3CYb3,OnlyDiagonal) 
 adjTYw3TYw3TpYw3CYw3 = Matmul2(adjTYw3,TYw3TpYw3CYw3,OnlyDiagonal) 
 adjTYx3TYx3TpYx3CYx3 = Matmul2(adjTYx3,TYx3TpYx3CYx3,OnlyDiagonal) 
 CYb3ml2TpYb3CYb3 = Matmul2(Conjg(Yb3),ml2TpYb3CYb3,OnlyDiagonal) 
 CYb3ml2TpYw3CYw3 = Matmul2(Conjg(Yb3),ml2TpYw3CYw3,OnlyDiagonal) 
 CYb3TpYb3mHb32CYb3 = Matmul2(Conjg(Yb3),TpYb3mHb32CYb3,OnlyDiagonal) 
 CYb3TpYb3CYb3ml2 = Matmul2(Conjg(Yb3),TpYb3CYb3ml2,OnlyDiagonal) 
 CYb3TpYb3CYb3TpYb3 = Matmul2(Conjg(Yb3),TpYb3CYb3TpYb3,OnlyDiagonal) 
Forall(i2=1:3)  CYb3TpYb3CYb3TpYb3(i2,i2) =  Real(CYb3TpYb3CYb3TpYb3(i2,i2),dp) 
 CYb3TpYb3CYb3TpTYb3 = Matmul2(Conjg(Yb3),TpYb3CYb3TpTYb3,OnlyDiagonal) 
 CYb3TpYeCYeTpYb3 = Matmul2(Conjg(Yb3),TpYeCYeTpYb3,OnlyDiagonal) 
Forall(i2=1:3)  CYb3TpYeCYeTpYb3(i2,i2) =  Real(CYb3TpYeCYeTpYb3(i2,i2),dp) 
 CYb3TpYeCYeTpTYb3 = Matmul2(Conjg(Yb3),TpYeCYeTpTYb3,OnlyDiagonal) 
 CYb3TpYw3mHw32CYw3 = Matmul2(Conjg(Yb3),TpYw3mHw32CYw3,OnlyDiagonal) 
 CYb3TpYw3CYw3ml2 = Matmul2(Conjg(Yb3),TpYw3CYw3ml2,OnlyDiagonal) 
 CYb3TpYw3CYw3TpYb3 = Matmul2(Conjg(Yb3),TpYw3CYw3TpYb3,OnlyDiagonal) 
Forall(i2=1:3)  CYb3TpYw3CYw3TpYb3(i2,i2) =  Real(CYb3TpYw3CYw3TpYb3(i2,i2),dp) 
 CYb3TpYw3CYw3TpTYb3 = Matmul2(Conjg(Yb3),TpYw3CYw3TpTYb3,OnlyDiagonal) 
 CYb3TpTYb3CYb3TpYb3 = Matmul2(Conjg(Yb3),TpTYb3CYb3TpYb3,OnlyDiagonal) 
 CYb3TpTYeCYeTpYb3 = Matmul2(Conjg(Yb3),TpTYeCYeTpYb3,OnlyDiagonal) 
 CYb3TpTYw3CYw3TpYb3 = Matmul2(Conjg(Yb3),TpTYw3CYw3TpYb3,OnlyDiagonal) 
 CYdmq2TpYdadjYx3 = Matmul2(Conjg(Yd),mq2TpYdadjYx3,OnlyDiagonal) 
 CYdmq2TpYdCYd = Matmul2(Conjg(Yd),mq2TpYdCYd,OnlyDiagonal) 
 CYdTpYdmd2adjYx3 = Matmul2(Conjg(Yd),TpYdmd2adjYx3,OnlyDiagonal) 
 CYdTpYdmd2CYd = Matmul2(Conjg(Yd),TpYdmd2CYd,OnlyDiagonal) 
 CYdTpYdadjYx3mHxb32 = Matmul2(Conjg(Yd),TpYdadjYx3mHxb32,OnlyDiagonal) 
 CYdTpYdadjYx3Yx3 = Matmul2(Conjg(Yd),TpYdadjYx3Yx3,OnlyDiagonal) 
 CYdTpYdadjYx3TYx3 = Matmul2(Conjg(Yd),TpYdadjYx3TYx3,OnlyDiagonal) 
 CYdTpYdCYdmq2 = Matmul2(Conjg(Yd),TpYdCYdmq2,OnlyDiagonal) 
 CYdTpYdCYdTpYd = Matmul2(Conjg(Yd),TpYdCYdTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYdCYdTpYd(i2,i2) =  Real(CYdTpYdCYdTpYd(i2,i2),dp) 
 CYdTpYdCYdTpTYd = Matmul2(Conjg(Yd),TpYdCYdTpTYd,OnlyDiagonal) 
 CYdTpYuCYuTpYd = Matmul2(Conjg(Yd),TpYuCYuTpYd,OnlyDiagonal) 
Forall(i2=1:3)  CYdTpYuCYuTpYd(i2,i2) =  Real(CYdTpYuCYuTpYd(i2,i2),dp) 
 CYdTpYuCYuTpTYd = Matmul2(Conjg(Yd),TpYuCYuTpTYd,OnlyDiagonal) 
 CYdTpTYdCYdTpYd = Matmul2(Conjg(Yd),TpTYdCYdTpYd,OnlyDiagonal) 
 CYdTpTYuCYuTpYd = Matmul2(Conjg(Yd),TpTYuCYuTpYd,OnlyDiagonal) 
 CYeml2TpYeCYe = Matmul2(Conjg(Ye),ml2TpYeCYe,OnlyDiagonal) 
 CYeTpYeme2CYe = Matmul2(Conjg(Ye),TpYeme2CYe,OnlyDiagonal) 
 CYeTpYeCYeml2 = Matmul2(Conjg(Ye),TpYeCYeml2,OnlyDiagonal) 
 CYumq2TpYuCYu = Matmul2(Conjg(Yu),mq2TpYuCYu,OnlyDiagonal) 
 CYuTpYumu2CYu = Matmul2(Conjg(Yu),TpYumu2CYu,OnlyDiagonal) 
 CYuTpYuCYumq2 = Matmul2(Conjg(Yu),TpYuCYumq2,OnlyDiagonal) 
 CYw3ml2TpYb3CYb3 = Matmul2(Conjg(Yw3),ml2TpYb3CYb3,OnlyDiagonal) 
 CYw3ml2TpYw3CYw3 = Matmul2(Conjg(Yw3),ml2TpYw3CYw3,OnlyDiagonal) 
 CYw3TpYb3mHb32CYb3 = Matmul2(Conjg(Yw3),TpYb3mHb32CYb3,OnlyDiagonal) 
 CYw3TpYb3CYb3ml2 = Matmul2(Conjg(Yw3),TpYb3CYb3ml2,OnlyDiagonal) 
 CYw3TpYb3CYb3TpYw3 = Matmul2(Conjg(Yw3),TpYb3CYb3TpYw3,OnlyDiagonal) 
Forall(i2=1:3)  CYw3TpYb3CYb3TpYw3(i2,i2) =  Real(CYw3TpYb3CYb3TpYw3(i2,i2),dp) 
 CYw3TpYb3CYb3TpTYw3 = Matmul2(Conjg(Yw3),TpYb3CYb3TpTYw3,OnlyDiagonal) 
 CYw3TpYeCYeTpYw3 = Matmul2(Conjg(Yw3),TpYeCYeTpYw3,OnlyDiagonal) 
Forall(i2=1:3)  CYw3TpYeCYeTpYw3(i2,i2) =  Real(CYw3TpYeCYeTpYw3(i2,i2),dp) 
 CYw3TpYeCYeTpTYw3 = Matmul2(Conjg(Yw3),TpYeCYeTpTYw3,OnlyDiagonal) 
 CYw3TpYw3mHw32CYw3 = Matmul2(Conjg(Yw3),TpYw3mHw32CYw3,OnlyDiagonal) 
 CYw3TpYw3CYw3ml2 = Matmul2(Conjg(Yw3),TpYw3CYw3ml2,OnlyDiagonal) 
 CYw3TpYw3CYw3TpYw3 = Matmul2(Conjg(Yw3),TpYw3CYw3TpYw3,OnlyDiagonal) 
Forall(i2=1:3)  CYw3TpYw3CYw3TpYw3(i2,i2) =  Real(CYw3TpYw3CYw3TpYw3(i2,i2),dp) 
 CYw3TpYw3CYw3TpTYw3 = Matmul2(Conjg(Yw3),TpYw3CYw3TpTYw3,OnlyDiagonal) 
 CYw3TpTYb3CYb3TpYw3 = Matmul2(Conjg(Yw3),TpTYb3CYb3TpYw3,OnlyDiagonal) 
 CYw3TpTYeCYeTpYw3 = Matmul2(Conjg(Yw3),TpTYeCYeTpYw3,OnlyDiagonal) 
 CYw3TpTYw3CYw3TpYw3 = Matmul2(Conjg(Yw3),TpTYw3CYw3TpYw3,OnlyDiagonal) 
 CYx3md2YdadjYd = Matmul2(Conjg(Yx3),md2YdadjYd,OnlyDiagonal) 
 CYx3md2TpYx3CYx3 = Matmul2(Conjg(Yx3),md2TpYx3CYx3,OnlyDiagonal) 
 CYx3YdadjYdTpYx3 = Matmul2(Conjg(Yx3),YdadjYdTpYx3,OnlyDiagonal) 
Forall(i2=1:3)  CYx3YdadjYdTpYx3(i2,i2) =  Real(CYx3YdadjYdTpYx3(i2,i2),dp) 
 CYx3YdadjYdTpTYx3 = Matmul2(Conjg(Yx3),YdadjYdTpTYx3,OnlyDiagonal) 
 CYx3TYdadjYdTpYx3 = Matmul2(Conjg(Yx3),TYdadjYdTpYx3,OnlyDiagonal) 
 CYx3TpYx3mHxb32CYx3 = Matmul2(Conjg(Yx3),TpYx3mHxb32CYx3,OnlyDiagonal) 
 CYx3TpYx3CYx3md2 = Matmul2(Conjg(Yx3),TpYx3CYx3md2,OnlyDiagonal) 
 CYx3TpYx3CYx3Yd = Matmul2(Conjg(Yx3),TpYx3CYx3Yd,OnlyDiagonal) 
 CYx3TpYx3CYx3TYd = Matmul2(Conjg(Yx3),TpYx3CYx3TYd,OnlyDiagonal) 
 CYx3TpYx3CYx3TpYx3 = Matmul2(Conjg(Yx3),TpYx3CYx3TpYx3,OnlyDiagonal) 
Forall(i2=1:3)  CYx3TpYx3CYx3TpYx3(i2,i2) =  Real(CYx3TpYx3CYx3TpYx3(i2,i2),dp) 
 CYx3TpYx3CYx3TpTYx3 = Matmul2(Conjg(Yx3),TpYx3CYx3TpTYx3,OnlyDiagonal) 
 CYx3TpTYx3CYx3Yd = Matmul2(Conjg(Yx3),TpTYx3CYx3Yd,OnlyDiagonal) 
 CYx3TpTYx3CYx3TpYx3 = Matmul2(Conjg(Yx3),TpTYx3CYx3TpYx3,OnlyDiagonal) 
 CTYdTpYdadjYx3TYx3 = Matmul2(Conjg(TYd),TpYdadjYx3TYx3,OnlyDiagonal) 
 TYb3adjYb3Yb3adjTYb3 = Matmul2(TYb3,adjYb3Yb3adjTYb3,OnlyDiagonal) 
 TYb3adjYeYeadjTYb3 = Matmul2(TYb3,adjYeYeadjTYb3,OnlyDiagonal) 
 TYb3adjYw3Yw3adjTYb3 = Matmul2(TYb3,adjYw3Yw3adjTYb3,OnlyDiagonal) 
 TYb3TpYb3CTYb3adjYb3 = Matmul2(TYb3,TpYb3CTYb3adjYb3,OnlyDiagonal) 
 TYb3TpYeCTYeadjYb3 = Matmul2(TYb3,TpYeCTYeadjYb3,OnlyDiagonal) 
 TYb3TpYw3CTYw3adjYb3 = Matmul2(TYb3,TpYw3CTYw3adjYb3,OnlyDiagonal) 
 TYdadjYdYdadjTYd = Matmul2(TYd,adjYdYdadjTYd,OnlyDiagonal) 
 TYdadjYdadjYx3Yx3 = Matmul2(TYd,adjYdadjYx3Yx3,OnlyDiagonal) 
 TYdadjYdTpYx3CYx3 = Matmul2(TYd,adjYdTpYx3CYx3,OnlyDiagonal) 
 TYdadjYdTpYx3CTYx3 = Matmul2(TYd,adjYdTpYx3CTYx3,OnlyDiagonal) 
 TYdadjYuYuadjTYd = Matmul2(TYd,adjYuYuadjTYd,OnlyDiagonal) 
 TYdTpYdCTYdadjYd = Matmul2(TYd,TpYdCTYdadjYd,OnlyDiagonal) 
 TYdTpYuCTYuadjYd = Matmul2(TYd,TpYuCTYuadjYd,OnlyDiagonal) 
 TYeadjYb3Yb3adjTYe = Matmul2(TYe,adjYb3Yb3adjTYe,OnlyDiagonal) 
 TYeadjYeYeadjTYe = Matmul2(TYe,adjYeYeadjTYe,OnlyDiagonal) 
 TYeadjYw3Yw3adjTYe = Matmul2(TYe,adjYw3Yw3adjTYe,OnlyDiagonal) 
 TYeTpYb3CTYb3adjYe = Matmul2(TYe,TpYb3CTYb3adjYe,OnlyDiagonal) 
 TYeTpYeCTYeadjYe = Matmul2(TYe,TpYeCTYeadjYe,OnlyDiagonal) 
 TYeTpYw3CTYw3adjYe = Matmul2(TYe,TpYw3CTYw3adjYe,OnlyDiagonal) 
 TYuadjYdYdadjTYu = Matmul2(TYu,adjYdYdadjTYu,OnlyDiagonal) 
 TYuadjYuYuadjTYu = Matmul2(TYu,adjYuYuadjTYu,OnlyDiagonal) 
 TYuTpYdCTYdadjYu = Matmul2(TYu,TpYdCTYdadjYu,OnlyDiagonal) 
 TYuTpYuCTYuadjYu = Matmul2(TYu,TpYuCTYuadjYu,OnlyDiagonal) 
 TYw3adjYb3Yb3adjTYw3 = Matmul2(TYw3,adjYb3Yb3adjTYw3,OnlyDiagonal) 
 TYw3adjYeYeadjTYw3 = Matmul2(TYw3,adjYeYeadjTYw3,OnlyDiagonal) 
 TYw3adjYw3Yw3adjTYw3 = Matmul2(TYw3,adjYw3Yw3adjTYw3,OnlyDiagonal) 
 TYw3TpYb3CTYb3adjYw3 = Matmul2(TYw3,TpYb3CTYb3adjYw3,OnlyDiagonal) 
 TYw3TpYeCTYeadjYw3 = Matmul2(TYw3,TpYeCTYeadjYw3,OnlyDiagonal) 
 TYw3TpYw3CTYw3adjYw3 = Matmul2(TYw3,TpYw3CTYw3adjYw3,OnlyDiagonal) 
 TYx3YdadjTYdadjYx3 = Matmul2(TYx3,YdadjTYdadjYx3,OnlyDiagonal) 
 TYx3adjYx3Yx3adjTYx3 = Matmul2(TYx3,adjYx3Yx3adjTYx3,OnlyDiagonal) 
 TYx3CYdTpYdadjTYx3 = Matmul2(TYx3,CYdTpYdadjTYx3,OnlyDiagonal) 
 TYx3TpYx3CTYx3adjYx3 = Matmul2(TYx3,TpYx3CTYx3adjYx3,OnlyDiagonal) 
 TpYb3CYb3TpYb3CYb3 = Matmul2(Transpose(Yb3),CYb3TpYb3CYb3,OnlyDiagonal) 
Forall(i2=1:3)  TpYb3CYb3TpYb3CYb3(i2,i2) =  Real(TpYb3CYb3TpYb3CYb3(i2,i2),dp) 
 TpYb3CYb3TpYw3CYw3 = Matmul2(Transpose(Yb3),CYb3TpYw3CYw3,OnlyDiagonal) 
 TpYb3CYb3TpTYb3CTYb3 = Matmul2(Transpose(Yb3),CYb3TpTYb3CTYb3,OnlyDiagonal) 
 TpYb3CYb3TpTYw3CTYw3 = Matmul2(Transpose(Yb3),CYb3TpTYw3CTYw3,OnlyDiagonal) 
 TpYb3CTYb3TpTYb3CYb3 = Matmul2(Transpose(Yb3),CTYb3TpTYb3CYb3,OnlyDiagonal) 
 TpYb3CTYb3TpTYw3CYw3 = Matmul2(Transpose(Yb3),CTYb3TpTYw3CYw3,OnlyDiagonal) 
 TpYdadjYdTpTYx3CTYx3 = Matmul2(Transpose(Yd),adjYdTpTYx3CTYx3,OnlyDiagonal) 
 TpYdadjYx3mHxb32Yx3 = Matmul2(Transpose(Yd),adjYx3mHxb32Yx3,OnlyDiagonal) 
 TpYdadjYx3Yx3CYd = Matmul2(Transpose(Yd),adjYx3Yx3CYd,OnlyDiagonal) 
Forall(i2=1:3)  TpYdadjYx3Yx3CYd(i2,i2) =  Real(TpYdadjYx3Yx3CYd(i2,i2),dp) 
 TpYdadjYx3TYx3CTYd = Matmul2(Transpose(Yd),adjYx3TYx3CTYd,OnlyDiagonal) 
 TpYdCYdTpYdCYd = Matmul2(Transpose(Yd),CYdTpYdCYd,OnlyDiagonal) 
Forall(i2=1:3)  TpYdCYdTpYdCYd(i2,i2) =  Real(TpYdCYdTpYdCYd(i2,i2),dp) 
 TpYdCYdTpTYdCTYd = Matmul2(Transpose(Yd),CYdTpTYdCTYd,OnlyDiagonal) 
 TpYdCTYdTpTYdCYd = Matmul2(Transpose(Yd),CTYdTpTYdCYd,OnlyDiagonal) 
 TpYeCYeTpYeCYe = Matmul2(Transpose(Ye),CYeTpYeCYe,OnlyDiagonal) 
Forall(i2=1:3)  TpYeCYeTpYeCYe(i2,i2) =  Real(TpYeCYeTpYeCYe(i2,i2),dp) 
 TpYeCYeTpTYeCTYe = Matmul2(Transpose(Ye),CYeTpTYeCTYe,OnlyDiagonal) 
 TpYeCTYeTpTYeCYe = Matmul2(Transpose(Ye),CTYeTpTYeCYe,OnlyDiagonal) 
 TpYuCYuTpYuCYu = Matmul2(Transpose(Yu),CYuTpYuCYu,OnlyDiagonal) 
Forall(i2=1:3)  TpYuCYuTpYuCYu(i2,i2) =  Real(TpYuCYuTpYuCYu(i2,i2),dp) 
 TpYuCYuTpTYuCTYu = Matmul2(Transpose(Yu),CYuTpTYuCTYu,OnlyDiagonal) 
 TpYuCTYuTpTYuCYu = Matmul2(Transpose(Yu),CTYuTpTYuCYu,OnlyDiagonal) 
 TpYw3CYw3TpYb3CYb3 = Matmul2(Transpose(Yw3),CYw3TpYb3CYb3,OnlyDiagonal) 
 TpYw3CYw3TpYw3CYw3 = Matmul2(Transpose(Yw3),CYw3TpYw3CYw3,OnlyDiagonal) 
Forall(i2=1:3)  TpYw3CYw3TpYw3CYw3(i2,i2) =  Real(TpYw3CYw3TpYw3CYw3(i2,i2),dp) 
 TpYw3CYw3TpTYb3CTYb3 = Matmul2(Transpose(Yw3),CYw3TpTYb3CTYb3,OnlyDiagonal) 
 TpYw3CYw3TpTYw3CTYw3 = Matmul2(Transpose(Yw3),CYw3TpTYw3CTYw3,OnlyDiagonal) 
 TpYw3CTYw3TpTYb3CYb3 = Matmul2(Transpose(Yw3),CTYw3TpTYb3CYb3,OnlyDiagonal) 
 TpYw3CTYw3TpTYw3CYw3 = Matmul2(Transpose(Yw3),CTYw3TpTYw3CYw3,OnlyDiagonal) 
 TpYx3CYx3YdadjYd = Matmul2(Transpose(Yx3),CYx3YdadjYd,OnlyDiagonal) 
 TpYx3CYx3TYdadjYd = Matmul2(Transpose(Yx3),CYx3TYdadjYd,OnlyDiagonal) 
 TpYx3CYx3TpYx3CYx3 = Matmul2(Transpose(Yx3),CYx3TpYx3CYx3,OnlyDiagonal) 
Forall(i2=1:3)  TpYx3CYx3TpYx3CYx3(i2,i2) =  Real(TpYx3CYx3TpYx3CYx3(i2,i2),dp) 
 TpYx3CYx3TpTYx3CTYx3 = Matmul2(Transpose(Yx3),CYx3TpTYx3CTYx3,OnlyDiagonal) 
 TpYx3CTYx3TpTYx3CYx3 = Matmul2(Transpose(Yx3),CTYx3TpTYx3CYx3,OnlyDiagonal) 
 TpTYb3CYb3TpYb3CTYb3 = Matmul2(Transpose(TYb3),CYb3TpYb3CTYb3,OnlyDiagonal) 
 TpTYb3CYb3TpYw3CTYw3 = Matmul2(Transpose(TYb3),CYb3TpYw3CTYw3,OnlyDiagonal) 
 TpTYdadjYdTpYx3CTYx3 = Matmul2(Transpose(TYd),adjYdTpYx3CTYx3,OnlyDiagonal) 
 TpTYdadjYx3Yx3CTYd = Matmul2(Transpose(TYd),adjYx3Yx3CTYd,OnlyDiagonal) 
 TpTYdCYdTpYdCTYd = Matmul2(Transpose(TYd),CYdTpYdCTYd,OnlyDiagonal) 
 TpTYeCYeTpYeCTYe = Matmul2(Transpose(TYe),CYeTpYeCTYe,OnlyDiagonal) 
 TpTYuCYuTpYuCTYu = Matmul2(Transpose(TYu),CYuTpYuCTYu,OnlyDiagonal) 
 TpTYw3CYw3TpYb3CTYb3 = Matmul2(Transpose(TYw3),CYw3TpYb3CTYb3,OnlyDiagonal) 
 TpTYw3CYw3TpYw3CTYw3 = Matmul2(Transpose(TYw3),CYw3TpYw3CTYw3,OnlyDiagonal) 
 TpTYx3CYx3TpYx3CTYx3 = Matmul2(Transpose(TYx3),CYx3TpYx3CTYx3,OnlyDiagonal) 
 MBM3CYb3TpYb3CYb3TpYb3 = Matmul2(MBM3,CYb3TpYb3CYb3TpYb3,OnlyDiagonal) 
 MBM3CYb3TpYb3CYb3TpTYb3 = Matmul2(MBM3,CYb3TpYb3CYb3TpTYb3,OnlyDiagonal) 
 MBM3CYb3TpYeCYeTpYb3 = Matmul2(MBM3,CYb3TpYeCYeTpYb3,OnlyDiagonal) 
 MBM3CYb3TpYeCYeTpTYb3 = Matmul2(MBM3,CYb3TpYeCYeTpTYb3,OnlyDiagonal) 
 MBM3CYb3TpYw3CYw3TpYb3 = Matmul2(MBM3,CYb3TpYw3CYw3TpYb3,OnlyDiagonal) 
 MBM3CYb3TpYw3CYw3TpTYb3 = Matmul2(MBM3,CYb3TpYw3CYw3TpTYb3,OnlyDiagonal) 
 MBM3CYb3TpTYb3CYb3TpYb3 = Matmul2(MBM3,CYb3TpTYb3CYb3TpYb3,OnlyDiagonal) 
 MBM3CYb3TpTYeCYeTpYb3 = Matmul2(MBM3,CYb3TpTYeCYeTpYb3,OnlyDiagonal) 
 MBM3CYb3TpTYw3CYw3TpYb3 = Matmul2(MBM3,CYb3TpTYw3CYw3TpYb3,OnlyDiagonal) 
 md2YdadjYdYdadjYd = Matmul2(md2,YdadjYdYdadjYd,OnlyDiagonal) 
 md2YdadjYdTpYx3CYx3 = Matmul2(md2,YdadjYdTpYx3CYx3,OnlyDiagonal) 
 md2YdadjYuYuadjYd = Matmul2(md2,YdadjYuYuadjYd,OnlyDiagonal) 
 md2TpYx3CYx3YdadjYd = Matmul2(md2,TpYx3CYx3YdadjYd,OnlyDiagonal) 
 md2TpYx3CYx3TpYx3CYx3 = Matmul2(md2,TpYx3CYx3TpYx3CYx3,OnlyDiagonal) 
 me2YeadjYb3Yb3adjYe = Matmul2(me2,YeadjYb3Yb3adjYe,OnlyDiagonal) 
 me2YeadjYeYeadjYe = Matmul2(me2,YeadjYeYeadjYe,OnlyDiagonal) 
 me2YeadjYw3Yw3adjYe = Matmul2(me2,YeadjYw3Yw3adjYe,OnlyDiagonal) 
 mHb32Yb3adjYb3Yb3adjYb3 = Matmul2(mHb32,Yb3adjYb3Yb3adjYb3,OnlyDiagonal) 
 mHb32Yb3adjYeYeadjYb3 = Matmul2(mHb32,Yb3adjYeYeadjYb3,OnlyDiagonal) 
 mHb32Yb3adjYw3Yw3adjYb3 = Matmul2(mHb32,Yb3adjYw3Yw3adjYb3,OnlyDiagonal) 
 mHw32Yw3adjYb3Yb3adjYw3 = Matmul2(mHw32,Yw3adjYb3Yb3adjYw3,OnlyDiagonal) 
 mHw32Yw3adjYeYeadjYw3 = Matmul2(mHw32,Yw3adjYeYeadjYw3,OnlyDiagonal) 
 mHw32Yw3adjYw3Yw3adjYw3 = Matmul2(mHw32,Yw3adjYw3Yw3adjYw3,OnlyDiagonal) 
 mHxb32Yx3adjYx3Yx3adjYx3 = Matmul2(mHxb32,Yx3adjYx3Yx3adjYx3,OnlyDiagonal) 
 mHxb32Yx3CYdTpYdadjYx3 = Matmul2(mHxb32,Yx3CYdTpYdadjYx3,OnlyDiagonal) 
 ml2TpYb3CYb3TpYb3CYb3 = Matmul2(ml2,TpYb3CYb3TpYb3CYb3,OnlyDiagonal) 
 ml2TpYb3CYb3TpYw3CYw3 = Matmul2(ml2,TpYb3CYb3TpYw3CYw3,OnlyDiagonal) 
 ml2TpYeCYeTpYeCYe = Matmul2(ml2,TpYeCYeTpYeCYe,OnlyDiagonal) 
 ml2TpYw3CYw3TpYb3CYb3 = Matmul2(ml2,TpYw3CYw3TpYb3CYb3,OnlyDiagonal) 
 ml2TpYw3CYw3TpYw3CYw3 = Matmul2(ml2,TpYw3CYw3TpYw3CYw3,OnlyDiagonal) 
 mq2adjYdTpYx3CYx3Yd = Matmul2(mq2,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 mq2TpYdCYdTpYdCYd = Matmul2(mq2,TpYdCYdTpYdCYd,OnlyDiagonal) 
 mq2TpYuCYuTpYuCYu = Matmul2(mq2,TpYuCYuTpYuCYu,OnlyDiagonal) 
 mu2YuadjYdYdadjYu = Matmul2(mu2,YuadjYdYdadjYu,OnlyDiagonal) 
 mu2YuadjYuYuadjYu = Matmul2(mu2,YuadjYuYuadjYu,OnlyDiagonal) 
 MWM3CYw3TpYb3CYb3TpYw3 = Matmul2(MWM3,CYw3TpYb3CYb3TpYw3,OnlyDiagonal) 
 MWM3CYw3TpYb3CYb3TpTYw3 = Matmul2(MWM3,CYw3TpYb3CYb3TpTYw3,OnlyDiagonal) 
 MWM3CYw3TpYeCYeTpYw3 = Matmul2(MWM3,CYw3TpYeCYeTpYw3,OnlyDiagonal) 
 MWM3CYw3TpYeCYeTpTYw3 = Matmul2(MWM3,CYw3TpYeCYeTpTYw3,OnlyDiagonal) 
 MWM3CYw3TpYw3CYw3TpYw3 = Matmul2(MWM3,CYw3TpYw3CYw3TpYw3,OnlyDiagonal) 
 MWM3CYw3TpYw3CYw3TpTYw3 = Matmul2(MWM3,CYw3TpYw3CYw3TpTYw3,OnlyDiagonal) 
 MWM3CYw3TpTYb3CYb3TpYw3 = Matmul2(MWM3,CYw3TpTYb3CYb3TpYw3,OnlyDiagonal) 
 MWM3CYw3TpTYeCYeTpYw3 = Matmul2(MWM3,CYw3TpTYeCYeTpYw3,OnlyDiagonal) 
 MWM3CYw3TpTYw3CYw3TpYw3 = Matmul2(MWM3,CYw3TpTYw3CYw3TpYw3,OnlyDiagonal) 
 MXM3CYx3YdadjYdTpYx3 = Matmul2(MXM3,CYx3YdadjYdTpYx3,OnlyDiagonal) 
 MXM3CYx3YdadjYdTpTYx3 = Matmul2(MXM3,CYx3YdadjYdTpTYx3,OnlyDiagonal) 
 MXM3CYx3TYdadjYdTpYx3 = Matmul2(MXM3,CYx3TYdadjYdTpYx3,OnlyDiagonal) 
 MXM3CYx3TpYx3CYx3TpYx3 = Matmul2(MXM3,CYx3TpYx3CYx3TpYx3,OnlyDiagonal) 
 MXM3CYx3TpYx3CYx3TpTYx3 = Matmul2(MXM3,CYx3TpYx3CYx3TpTYx3,OnlyDiagonal) 
 MXM3CYx3TpTYx3CYx3TpYx3 = Matmul2(MXM3,CYx3TpTYx3CYx3TpYx3,OnlyDiagonal) 
 Yb3ml2adjYb3Yb3adjYb3 = Matmul2(Yb3,ml2adjYb3Yb3adjYb3,OnlyDiagonal) 
 Yb3ml2adjYeYeadjYb3 = Matmul2(Yb3,ml2adjYeYeadjYb3,OnlyDiagonal) 
 Yb3ml2adjYw3Yw3adjYb3 = Matmul2(Yb3,ml2adjYw3Yw3adjYb3,OnlyDiagonal) 
 Yb3adjYb3mHb32Yb3adjYb3 = Matmul2(Yb3,adjYb3mHb32Yb3adjYb3,OnlyDiagonal) 
 Yb3adjYb3Yb3ml2adjYb3 = Matmul2(Yb3,adjYb3Yb3ml2adjYb3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3MBM3 = Matmul2(Yb3,adjYb3Yb3adjYb3MBM3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3mHb32 = Matmul2(Yb3,adjYb3Yb3adjYb3mHb32,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3Yb3 = Matmul2(Yb3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3BMBM3 = Matmul2(Yb3,adjYb3Yb3adjYb3BMBM3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYb3TYb3 = Matmul2(Yb3,adjYb3Yb3adjYb3TYb3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYw3Yw3 = Matmul2(Yb3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 Yb3adjYb3Yb3adjYw3TYw3 = Matmul2(Yb3,adjYb3Yb3adjYw3TYw3,OnlyDiagonal) 
 Yb3adjYb3TYb3adjYb3MBM3 = Matmul2(Yb3,adjYb3TYb3adjYb3MBM3,OnlyDiagonal) 
 Yb3adjYb3TYb3adjYb3Yb3 = Matmul2(Yb3,adjYb3TYb3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYb3TYb3adjYw3Yw3 = Matmul2(Yb3,adjYb3TYb3adjYw3Yw3,OnlyDiagonal) 
 Yb3adjYeme2YeadjYb3 = Matmul2(Yb3,adjYeme2YeadjYb3,OnlyDiagonal) 
 Yb3adjYeYeml2adjYb3 = Matmul2(Yb3,adjYeYeml2adjYb3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3MBM3 = Matmul2(Yb3,adjYeYeadjYb3MBM3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3mHb32 = Matmul2(Yb3,adjYeYeadjYb3mHb32,OnlyDiagonal) 
 Yb3adjYeYeadjYb3Yb3 = Matmul2(Yb3,adjYeYeadjYb3Yb3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3BMBM3 = Matmul2(Yb3,adjYeYeadjYb3BMBM3,OnlyDiagonal) 
 Yb3adjYeYeadjYb3TYb3 = Matmul2(Yb3,adjYeYeadjYb3TYb3,OnlyDiagonal) 
 Yb3adjYeYeadjYeYe = Matmul2(Yb3,adjYeYeadjYeYe,OnlyDiagonal) 
 Yb3adjYeYeadjYeTYe = Matmul2(Yb3,adjYeYeadjYeTYe,OnlyDiagonal) 
 Yb3adjYeYeadjYw3TYw3 = Matmul2(Yb3,adjYeYeadjYw3TYw3,OnlyDiagonal) 
 Yb3adjYeTYeadjYb3MBM3 = Matmul2(Yb3,adjYeTYeadjYb3MBM3,OnlyDiagonal) 
 Yb3adjYeTYeadjYb3Yb3 = Matmul2(Yb3,adjYeTYeadjYb3Yb3,OnlyDiagonal) 
 Yb3adjYeTYeadjYeYe = Matmul2(Yb3,adjYeTYeadjYeYe,OnlyDiagonal) 
 Yb3adjYeTYeadjYw3Yw3 = Matmul2(Yb3,adjYeTYeadjYw3Yw3,OnlyDiagonal) 
 Yb3adjYw3mHw32Yw3adjYb3 = Matmul2(Yb3,adjYw3mHw32Yw3adjYb3,OnlyDiagonal) 
 Yb3adjYw3Yw3ml2adjYb3 = Matmul2(Yb3,adjYw3Yw3ml2adjYb3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3MBM3 = Matmul2(Yb3,adjYw3Yw3adjYb3MBM3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3mHb32 = Matmul2(Yb3,adjYw3Yw3adjYb3mHb32,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3Yb3 = Matmul2(Yb3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3BMBM3 = Matmul2(Yb3,adjYw3Yw3adjYb3BMBM3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYb3TYb3 = Matmul2(Yb3,adjYw3Yw3adjYb3TYb3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYw3Yw3 = Matmul2(Yb3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 Yb3adjYw3Yw3adjYw3TYw3 = Matmul2(Yb3,adjYw3Yw3adjYw3TYw3,OnlyDiagonal) 
 Yb3adjYw3TYw3adjYb3MBM3 = Matmul2(Yb3,adjYw3TYw3adjYb3MBM3,OnlyDiagonal) 
 Yb3adjYw3TYw3adjYb3Yb3 = Matmul2(Yb3,adjYw3TYw3adjYb3Yb3,OnlyDiagonal) 
 Yb3adjYw3TYw3adjYw3Yw3 = Matmul2(Yb3,adjYw3TYw3adjYw3Yw3,OnlyDiagonal) 
 Ydmq2adjYdYdadjYd = Matmul2(Yd,mq2adjYdYdadjYd,OnlyDiagonal) 
 Ydmq2adjYdTpYx3CYx3 = Matmul2(Yd,mq2adjYdTpYx3CYx3,OnlyDiagonal) 
 Ydmq2adjYuYuadjYd = Matmul2(Yd,mq2adjYuYuadjYd,OnlyDiagonal) 
 Ydmq2adjYx3Yx3CYd = Matmul2(Yd,mq2adjYx3Yx3CYd,OnlyDiagonal) 
 YdadjYdmd2YdadjYd = Matmul2(Yd,adjYdmd2YdadjYd,OnlyDiagonal) 
 YdadjYdmd2TpYx3CYx3 = Matmul2(Yd,adjYdmd2TpYx3CYx3,OnlyDiagonal) 
 YdadjYdYdmq2adjYd = Matmul2(Yd,adjYdYdmq2adjYd,OnlyDiagonal) 
 YdadjYdYdadjYdmd2 = Matmul2(Yd,adjYdYdadjYdmd2,OnlyDiagonal) 
 YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal) 
 YdadjYdYdadjYdTYd = Matmul2(Yd,adjYdYdadjYdTYd,OnlyDiagonal) 
 YdadjYdTYdadjYdYd = Matmul2(Yd,adjYdTYdadjYdYd,OnlyDiagonal) 
 YdadjYdTpYx3CYx3Yd = Matmul2(Yd,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 YdadjYdTpYx3CYx3TYd = Matmul2(Yd,adjYdTpYx3CYx3TYd,OnlyDiagonal) 
 YdadjYdTpTYx3CYx3Yd = Matmul2(Yd,adjYdTpTYx3CYx3Yd,OnlyDiagonal) 
 YdadjYumu2YuadjYd = Matmul2(Yd,adjYumu2YuadjYd,OnlyDiagonal) 
 YdadjYuYumq2adjYd = Matmul2(Yd,adjYuYumq2adjYd,OnlyDiagonal) 
 YdadjYuYuadjYdmd2 = Matmul2(Yd,adjYuYuadjYdmd2,OnlyDiagonal) 
 YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal) 
 YdadjYuYuadjYdTYd = Matmul2(Yd,adjYuYuadjYdTYd,OnlyDiagonal) 
 YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal) 
 YdadjYuYuadjYuTYu = Matmul2(Yd,adjYuYuadjYuTYu,OnlyDiagonal) 
 YdadjYuTYuadjYdYd = Matmul2(Yd,adjYuTYuadjYdYd,OnlyDiagonal) 
 YdadjYuTYuadjYuYu = Matmul2(Yd,adjYuTYuadjYuYu,OnlyDiagonal) 
 Yeml2adjYb3Yb3adjYe = Matmul2(Ye,ml2adjYb3Yb3adjYe,OnlyDiagonal) 
 Yeml2adjYeYeadjYe = Matmul2(Ye,ml2adjYeYeadjYe,OnlyDiagonal) 
 Yeml2adjYw3Yw3adjYe = Matmul2(Ye,ml2adjYw3Yw3adjYe,OnlyDiagonal) 
 YeadjYb3mHb32Yb3adjYe = Matmul2(Ye,adjYb3mHb32Yb3adjYe,OnlyDiagonal) 
 YeadjYb3Yb3ml2adjYe = Matmul2(Ye,adjYb3Yb3ml2adjYe,OnlyDiagonal) 
 YeadjYb3Yb3adjYb3Yb3 = Matmul2(Ye,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 YeadjYb3Yb3adjYb3TYb3 = Matmul2(Ye,adjYb3Yb3adjYb3TYb3,OnlyDiagonal) 
 YeadjYb3Yb3adjYeme2 = Matmul2(Ye,adjYb3Yb3adjYeme2,OnlyDiagonal) 
 YeadjYb3Yb3adjYeYe = Matmul2(Ye,adjYb3Yb3adjYeYe,OnlyDiagonal) 
 YeadjYb3Yb3adjYeTYe = Matmul2(Ye,adjYb3Yb3adjYeTYe,OnlyDiagonal) 
 YeadjYb3Yb3adjYw3Yw3 = Matmul2(Ye,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 YeadjYb3Yb3adjYw3TYw3 = Matmul2(Ye,adjYb3Yb3adjYw3TYw3,OnlyDiagonal) 
 YeadjYb3TYb3adjYb3Yb3 = Matmul2(Ye,adjYb3TYb3adjYb3Yb3,OnlyDiagonal) 
 YeadjYb3TYb3adjYeYe = Matmul2(Ye,adjYb3TYb3adjYeYe,OnlyDiagonal) 
 YeadjYb3TYb3adjYw3Yw3 = Matmul2(Ye,adjYb3TYb3adjYw3Yw3,OnlyDiagonal) 
 YeadjYeme2YeadjYe = Matmul2(Ye,adjYeme2YeadjYe,OnlyDiagonal) 
 YeadjYeYeml2adjYe = Matmul2(Ye,adjYeYeml2adjYe,OnlyDiagonal) 
 YeadjYeYeadjYeme2 = Matmul2(Ye,adjYeYeadjYeme2,OnlyDiagonal) 
 YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal) 
 YeadjYeYeadjYeTYe = Matmul2(Ye,adjYeYeadjYeTYe,OnlyDiagonal) 
 YeadjYeTYeadjYeYe = Matmul2(Ye,adjYeTYeadjYeYe,OnlyDiagonal) 
 YeadjYw3mHw32Yw3adjYe = Matmul2(Ye,adjYw3mHw32Yw3adjYe,OnlyDiagonal) 
 YeadjYw3Yw3ml2adjYe = Matmul2(Ye,adjYw3Yw3ml2adjYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYb3Yb3 = Matmul2(Ye,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 YeadjYw3Yw3adjYb3TYb3 = Matmul2(Ye,adjYw3Yw3adjYb3TYb3,OnlyDiagonal) 
 YeadjYw3Yw3adjYeme2 = Matmul2(Ye,adjYw3Yw3adjYeme2,OnlyDiagonal) 
 YeadjYw3Yw3adjYeYe = Matmul2(Ye,adjYw3Yw3adjYeYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYeTYe = Matmul2(Ye,adjYw3Yw3adjYeTYe,OnlyDiagonal) 
 YeadjYw3Yw3adjYw3Yw3 = Matmul2(Ye,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 YeadjYw3Yw3adjYw3TYw3 = Matmul2(Ye,adjYw3Yw3adjYw3TYw3,OnlyDiagonal) 
 YeadjYw3TYw3adjYb3Yb3 = Matmul2(Ye,adjYw3TYw3adjYb3Yb3,OnlyDiagonal) 
 YeadjYw3TYw3adjYeYe = Matmul2(Ye,adjYw3TYw3adjYeYe,OnlyDiagonal) 
 YeadjYw3TYw3adjYw3Yw3 = Matmul2(Ye,adjYw3TYw3adjYw3Yw3,OnlyDiagonal) 
 Yumq2adjYdYdadjYu = Matmul2(Yu,mq2adjYdYdadjYu,OnlyDiagonal) 
 Yumq2adjYuYuadjYu = Matmul2(Yu,mq2adjYuYuadjYu,OnlyDiagonal) 
 YuadjYdmd2YdadjYu = Matmul2(Yu,adjYdmd2YdadjYu,OnlyDiagonal) 
 YuadjYdYdmq2adjYu = Matmul2(Yu,adjYdYdmq2adjYu,OnlyDiagonal) 
 YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal) 
 YuadjYdYdadjYdTYd = Matmul2(Yu,adjYdYdadjYdTYd,OnlyDiagonal) 
 YuadjYdYdadjYumu2 = Matmul2(Yu,adjYdYdadjYumu2,OnlyDiagonal) 
 YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal) 
 YuadjYdYdadjYuTYu = Matmul2(Yu,adjYdYdadjYuTYu,OnlyDiagonal) 
 YuadjYdTYdadjYdYd = Matmul2(Yu,adjYdTYdadjYdYd,OnlyDiagonal) 
 YuadjYdTYdadjYuYu = Matmul2(Yu,adjYdTYdadjYuYu,OnlyDiagonal) 
 YuadjYdTpYx3CYx3Yd = Matmul2(Yu,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 YuadjYdTpYx3CYx3TYd = Matmul2(Yu,adjYdTpYx3CYx3TYd,OnlyDiagonal) 
 YuadjYdTpTYx3CYx3Yd = Matmul2(Yu,adjYdTpTYx3CYx3Yd,OnlyDiagonal) 
 YuadjYumu2YuadjYu = Matmul2(Yu,adjYumu2YuadjYu,OnlyDiagonal) 
 YuadjYuYumq2adjYu = Matmul2(Yu,adjYuYumq2adjYu,OnlyDiagonal) 
 YuadjYuYuadjYumu2 = Matmul2(Yu,adjYuYuadjYumu2,OnlyDiagonal) 
 YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal) 
 YuadjYuYuadjYuTYu = Matmul2(Yu,adjYuYuadjYuTYu,OnlyDiagonal) 
 YuadjYuTYuadjYuYu = Matmul2(Yu,adjYuTYuadjYuYu,OnlyDiagonal) 
 Yw3ml2adjYb3Yb3adjYw3 = Matmul2(Yw3,ml2adjYb3Yb3adjYw3,OnlyDiagonal) 
 Yw3ml2adjYeYeadjYw3 = Matmul2(Yw3,ml2adjYeYeadjYw3,OnlyDiagonal) 
 Yw3ml2adjYw3Yw3adjYw3 = Matmul2(Yw3,ml2adjYw3Yw3adjYw3,OnlyDiagonal) 
 Yw3adjYb3mHb32Yb3adjYw3 = Matmul2(Yw3,adjYb3mHb32Yb3adjYw3,OnlyDiagonal) 
 Yw3adjYb3Yb3ml2adjYw3 = Matmul2(Yw3,adjYb3Yb3ml2adjYw3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYb3Yb3 = Matmul2(Yw3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYb3TYb3 = Matmul2(Yw3,adjYb3Yb3adjYb3TYb3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3mHw32 = Matmul2(Yw3,adjYb3Yb3adjYw3mHw32,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3MWM3 = Matmul2(Yw3,adjYb3Yb3adjYw3MWM3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3Yw3 = Matmul2(Yw3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3BMWM3 = Matmul2(Yw3,adjYb3Yb3adjYw3BMWM3,OnlyDiagonal) 
 Yw3adjYb3Yb3adjYw3TYw3 = Matmul2(Yw3,adjYb3Yb3adjYw3TYw3,OnlyDiagonal) 
 Yw3adjYb3TYb3adjYb3Yb3 = Matmul2(Yw3,adjYb3TYb3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYb3TYb3adjYw3MWM3 = Matmul2(Yw3,adjYb3TYb3adjYw3MWM3,OnlyDiagonal) 
 Yw3adjYb3TYb3adjYw3Yw3 = Matmul2(Yw3,adjYb3TYb3adjYw3Yw3,OnlyDiagonal) 
 Yw3adjYeme2YeadjYw3 = Matmul2(Yw3,adjYeme2YeadjYw3,OnlyDiagonal) 
 Yw3adjYeYeml2adjYw3 = Matmul2(Yw3,adjYeYeml2adjYw3,OnlyDiagonal) 
 Yw3adjYeYeadjYb3TYb3 = Matmul2(Yw3,adjYeYeadjYb3TYb3,OnlyDiagonal) 
 Yw3adjYeYeadjYeYe = Matmul2(Yw3,adjYeYeadjYeYe,OnlyDiagonal) 
 Yw3adjYeYeadjYeTYe = Matmul2(Yw3,adjYeYeadjYeTYe,OnlyDiagonal) 
 Yw3adjYeYeadjYw3mHw32 = Matmul2(Yw3,adjYeYeadjYw3mHw32,OnlyDiagonal) 
 Yw3adjYeYeadjYw3MWM3 = Matmul2(Yw3,adjYeYeadjYw3MWM3,OnlyDiagonal) 
 Yw3adjYeYeadjYw3Yw3 = Matmul2(Yw3,adjYeYeadjYw3Yw3,OnlyDiagonal) 
 Yw3adjYeYeadjYw3BMWM3 = Matmul2(Yw3,adjYeYeadjYw3BMWM3,OnlyDiagonal) 
 Yw3adjYeYeadjYw3TYw3 = Matmul2(Yw3,adjYeYeadjYw3TYw3,OnlyDiagonal) 
 Yw3adjYeTYeadjYb3Yb3 = Matmul2(Yw3,adjYeTYeadjYb3Yb3,OnlyDiagonal) 
 Yw3adjYeTYeadjYeYe = Matmul2(Yw3,adjYeTYeadjYeYe,OnlyDiagonal) 
 Yw3adjYeTYeadjYw3MWM3 = Matmul2(Yw3,adjYeTYeadjYw3MWM3,OnlyDiagonal) 
 Yw3adjYeTYeadjYw3Yw3 = Matmul2(Yw3,adjYeTYeadjYw3Yw3,OnlyDiagonal) 
 Yw3adjYw3mHw32Yw3adjYw3 = Matmul2(Yw3,adjYw3mHw32Yw3adjYw3,OnlyDiagonal) 
 Yw3adjYw3Yw3ml2adjYw3 = Matmul2(Yw3,adjYw3Yw3ml2adjYw3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYb3Yb3 = Matmul2(Yw3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYb3TYb3 = Matmul2(Yw3,adjYw3Yw3adjYb3TYb3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3mHw32 = Matmul2(Yw3,adjYw3Yw3adjYw3mHw32,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3MWM3 = Matmul2(Yw3,adjYw3Yw3adjYw3MWM3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3Yw3 = Matmul2(Yw3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3BMWM3 = Matmul2(Yw3,adjYw3Yw3adjYw3BMWM3,OnlyDiagonal) 
 Yw3adjYw3Yw3adjYw3TYw3 = Matmul2(Yw3,adjYw3Yw3adjYw3TYw3,OnlyDiagonal) 
 Yw3adjYw3TYw3adjYb3Yb3 = Matmul2(Yw3,adjYw3TYw3adjYb3Yb3,OnlyDiagonal) 
 Yw3adjYw3TYw3adjYw3MWM3 = Matmul2(Yw3,adjYw3TYw3adjYw3MWM3,OnlyDiagonal) 
 Yw3adjYw3TYw3adjYw3Yw3 = Matmul2(Yw3,adjYw3TYw3adjYw3Yw3,OnlyDiagonal) 
 Yx3md2adjYx3Yx3adjYx3 = Matmul2(Yx3,md2adjYx3Yx3adjYx3,OnlyDiagonal) 
 Yx3md2CYdTpYdadjYx3 = Matmul2(Yx3,md2CYdTpYdadjYx3,OnlyDiagonal) 
 Yx3adjYx3mHxb32Yx3adjYx3 = Matmul2(Yx3,adjYx3mHxb32Yx3adjYx3,OnlyDiagonal) 
 Yx3adjYx3Yx3md2adjYx3 = Matmul2(Yx3,adjYx3Yx3md2adjYx3,OnlyDiagonal) 
 Yx3adjYx3Yx3adjYx3mHxb32 = Matmul2(Yx3,adjYx3Yx3adjYx3mHxb32,OnlyDiagonal) 
 Yx3adjYx3Yx3adjYx3Yx3 = Matmul2(Yx3,adjYx3Yx3adjYx3Yx3,OnlyDiagonal) 
 Yx3adjYx3Yx3adjYx3TYx3 = Matmul2(Yx3,adjYx3Yx3adjYx3TYx3,OnlyDiagonal) 
 Yx3adjYx3TYx3adjYx3Yx3 = Matmul2(Yx3,adjYx3TYx3adjYx3Yx3,OnlyDiagonal) 
 Yx3CYdmq2TpYdadjYx3 = Matmul2(Yx3,CYdmq2TpYdadjYx3,OnlyDiagonal) 
 Yx3CYdTpYdmd2adjYx3 = Matmul2(Yx3,CYdTpYdmd2adjYx3,OnlyDiagonal) 
 Yx3CYdTpYdadjYx3mHxb32 = Matmul2(Yx3,CYdTpYdadjYx3mHxb32,OnlyDiagonal) 
 Yx3CYdTpYdadjYx3Yx3 = Matmul2(Yx3,CYdTpYdadjYx3Yx3,OnlyDiagonal) 
 Yx3CYdTpYdadjYx3TYx3 = Matmul2(Yx3,CYdTpYdadjYx3TYx3,OnlyDiagonal) 
 Yx3CYdTpYdCYdTpYd = Matmul2(Yx3,CYdTpYdCYdTpYd,OnlyDiagonal) 
 Yx3CYdTpYdCYdTpTYd = Matmul2(Yx3,CYdTpYdCYdTpTYd,OnlyDiagonal) 
 Yx3CYdTpYuCYuTpYd = Matmul2(Yx3,CYdTpYuCYuTpYd,OnlyDiagonal) 
 Yx3CYdTpYuCYuTpTYd = Matmul2(Yx3,CYdTpYuCYuTpTYd,OnlyDiagonal) 
 Yx3CYdTpTYdCYdTpYd = Matmul2(Yx3,CYdTpTYdCYdTpYd,OnlyDiagonal) 
 Yx3CYdTpTYuCYuTpYd = Matmul2(Yx3,CYdTpTYuCYuTpYd,OnlyDiagonal) 
 Yx3TYdadjYdadjYx3Yx3 = Matmul2(Yx3,TYdadjYdadjYx3Yx3,OnlyDiagonal) 
 BMBM3CYb3TpYb3CYb3TpYb3 = Matmul2(BMBM3,CYb3TpYb3CYb3TpYb3,OnlyDiagonal) 
 BMBM3CYb3TpYeCYeTpYb3 = Matmul2(BMBM3,CYb3TpYeCYeTpYb3,OnlyDiagonal) 
 BMBM3CYb3TpYw3CYw3TpYb3 = Matmul2(BMBM3,CYb3TpYw3CYw3TpYb3,OnlyDiagonal) 
 BMWM3CYw3TpYb3CYb3TpYw3 = Matmul2(BMWM3,CYw3TpYb3CYb3TpYw3,OnlyDiagonal) 
 BMWM3CYw3TpYeCYeTpYw3 = Matmul2(BMWM3,CYw3TpYeCYeTpYw3,OnlyDiagonal) 
 BMWM3CYw3TpYw3CYw3TpYw3 = Matmul2(BMWM3,CYw3TpYw3CYw3TpYw3,OnlyDiagonal) 
 BMXM3CYx3YdadjYdTpYx3 = Matmul2(BMXM3,CYx3YdadjYdTpYx3,OnlyDiagonal) 
 BMXM3CYx3TpYx3CYx3TpYx3 = Matmul2(BMXM3,CYx3TpYx3CYx3TpYx3,OnlyDiagonal) 
 CYdTpYdadjYx3mHxb32Yx3 = Matmul2(Conjg(Yd),TpYdadjYx3mHxb32Yx3,OnlyDiagonal) 
 TYb3adjYb3Yb3adjYb3MBM3 = Matmul2(TYb3,adjYb3Yb3adjYb3MBM3,OnlyDiagonal) 
 TYb3adjYb3Yb3adjYb3Yb3 = Matmul2(TYb3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 TYb3adjYb3Yb3adjYw3Yw3 = Matmul2(TYb3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 TYb3adjYeYeadjYb3MBM3 = Matmul2(TYb3,adjYeYeadjYb3MBM3,OnlyDiagonal) 
 TYb3adjYeYeadjYb3Yb3 = Matmul2(TYb3,adjYeYeadjYb3Yb3,OnlyDiagonal) 
 TYb3adjYeYeadjYeYe = Matmul2(TYb3,adjYeYeadjYeYe,OnlyDiagonal) 
 TYb3adjYeYeadjYw3Yw3 = Matmul2(TYb3,adjYeYeadjYw3Yw3,OnlyDiagonal) 
 TYb3adjYw3Yw3adjYb3MBM3 = Matmul2(TYb3,adjYw3Yw3adjYb3MBM3,OnlyDiagonal) 
 TYb3adjYw3Yw3adjYb3Yb3 = Matmul2(TYb3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 TYb3adjYw3Yw3adjYw3Yw3 = Matmul2(TYb3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 TYdadjYdYdadjYdYd = Matmul2(TYd,adjYdYdadjYdYd,OnlyDiagonal) 
 TYdadjYdTpYx3CYx3Yd = Matmul2(TYd,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 TYdadjYuYuadjYdYd = Matmul2(TYd,adjYuYuadjYdYd,OnlyDiagonal) 
 TYdadjYuYuadjYuYu = Matmul2(TYd,adjYuYuadjYuYu,OnlyDiagonal) 
 TYeadjYb3Yb3adjYb3Yb3 = Matmul2(TYe,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 TYeadjYb3Yb3adjYeYe = Matmul2(TYe,adjYb3Yb3adjYeYe,OnlyDiagonal) 
 TYeadjYb3Yb3adjYw3Yw3 = Matmul2(TYe,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 TYeadjYeYeadjYeYe = Matmul2(TYe,adjYeYeadjYeYe,OnlyDiagonal) 
 TYeadjYw3Yw3adjYb3Yb3 = Matmul2(TYe,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 TYeadjYw3Yw3adjYeYe = Matmul2(TYe,adjYw3Yw3adjYeYe,OnlyDiagonal) 
 TYeadjYw3Yw3adjYw3Yw3 = Matmul2(TYe,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 TYuadjYdYdadjYdYd = Matmul2(TYu,adjYdYdadjYdYd,OnlyDiagonal) 
 TYuadjYdYdadjYuYu = Matmul2(TYu,adjYdYdadjYuYu,OnlyDiagonal) 
 TYuadjYdTpYx3CYx3Yd = Matmul2(TYu,adjYdTpYx3CYx3Yd,OnlyDiagonal) 
 TYuadjYuYuadjYuYu = Matmul2(TYu,adjYuYuadjYuYu,OnlyDiagonal) 
 TYw3adjYb3Yb3adjYb3Yb3 = Matmul2(TYw3,adjYb3Yb3adjYb3Yb3,OnlyDiagonal) 
 TYw3adjYb3Yb3adjYw3MWM3 = Matmul2(TYw3,adjYb3Yb3adjYw3MWM3,OnlyDiagonal) 
 TYw3adjYb3Yb3adjYw3Yw3 = Matmul2(TYw3,adjYb3Yb3adjYw3Yw3,OnlyDiagonal) 
 TYw3adjYeYeadjYb3Yb3 = Matmul2(TYw3,adjYeYeadjYb3Yb3,OnlyDiagonal) 
 TYw3adjYeYeadjYeYe = Matmul2(TYw3,adjYeYeadjYeYe,OnlyDiagonal) 
 TYw3adjYeYeadjYw3MWM3 = Matmul2(TYw3,adjYeYeadjYw3MWM3,OnlyDiagonal) 
 TYw3adjYeYeadjYw3Yw3 = Matmul2(TYw3,adjYeYeadjYw3Yw3,OnlyDiagonal) 
 TYw3adjYw3Yw3adjYb3Yb3 = Matmul2(TYw3,adjYw3Yw3adjYb3Yb3,OnlyDiagonal) 
 TYw3adjYw3Yw3adjYw3MWM3 = Matmul2(TYw3,adjYw3Yw3adjYw3MWM3,OnlyDiagonal) 
 TYw3adjYw3Yw3adjYw3Yw3 = Matmul2(TYw3,adjYw3Yw3adjYw3Yw3,OnlyDiagonal) 
 TYx3adjYx3Yx3adjYx3Yx3 = Matmul2(TYx3,adjYx3Yx3adjYx3Yx3,OnlyDiagonal) 
 TYx3CYdTpYdadjYx3Yx3 = Matmul2(TYx3,CYdTpYdadjYx3Yx3,OnlyDiagonal) 
 TYx3CYdTpYdCYdTpYd = Matmul2(TYx3,CYdTpYdCYdTpYd,OnlyDiagonal) 
 TYx3CYdTpYuCYuTpYd = Matmul2(TYx3,CYdTpYuCYuTpYd,OnlyDiagonal) 
 TpYb3mHb32CYb3TpYb3CYb3 = Matmul2(Transpose(Yb3),mHb32CYb3TpYb3CYb3,OnlyDiagonal) 
 TpYb3mHb32CYb3TpYw3CYw3 = Matmul2(Transpose(Yb3),mHb32CYb3TpYw3CYw3,OnlyDiagonal) 
 TpYb3CYb3ml2TpYb3CYb3 = Matmul2(Transpose(Yb3),CYb3ml2TpYb3CYb3,OnlyDiagonal) 
 TpYb3CYb3ml2TpYw3CYw3 = Matmul2(Transpose(Yb3),CYb3ml2TpYw3CYw3,OnlyDiagonal) 
 TpYb3CYb3TpYb3mHb32CYb3 = Matmul2(Transpose(Yb3),CYb3TpYb3mHb32CYb3,OnlyDiagonal) 
 TpYb3CYb3TpYb3CYb3ml2 = Matmul2(Transpose(Yb3),CYb3TpYb3CYb3ml2,OnlyDiagonal) 
 TpYb3CYb3TpYw3mHw32CYw3 = Matmul2(Transpose(Yb3),CYb3TpYw3mHw32CYw3,OnlyDiagonal) 
 TpYb3CYb3TpYw3CYw3ml2 = Matmul2(Transpose(Yb3),CYb3TpYw3CYw3ml2,OnlyDiagonal) 
 TpYdmd2adjYx3Yx3CYd = Matmul2(Transpose(Yd),md2adjYx3Yx3CYd,OnlyDiagonal) 
 TpYdmd2CYdTpYdCYd = Matmul2(Transpose(Yd),md2CYdTpYdCYd,OnlyDiagonal) 
 TpYdadjYx3mHxb32Yx3CYd = Matmul2(Transpose(Yd),adjYx3mHxb32Yx3CYd,OnlyDiagonal) 
 TpYdadjYx3Yx3md2CYd = Matmul2(Transpose(Yd),adjYx3Yx3md2CYd,OnlyDiagonal) 
 TpYdadjYx3Yx3CYdmq2 = Matmul2(Transpose(Yd),adjYx3Yx3CYdmq2,OnlyDiagonal) 
 TpYdCYdmq2TpYdCYd = Matmul2(Transpose(Yd),CYdmq2TpYdCYd,OnlyDiagonal) 
 TpYdCYdTpYdmd2CYd = Matmul2(Transpose(Yd),CYdTpYdmd2CYd,OnlyDiagonal) 
 TpYdCYdTpYdCYdmq2 = Matmul2(Transpose(Yd),CYdTpYdCYdmq2,OnlyDiagonal) 
 TpYeme2CYeTpYeCYe = Matmul2(Transpose(Ye),me2CYeTpYeCYe,OnlyDiagonal) 
 TpYeCYeml2TpYeCYe = Matmul2(Transpose(Ye),CYeml2TpYeCYe,OnlyDiagonal) 
 TpYeCYeTpYeme2CYe = Matmul2(Transpose(Ye),CYeTpYeme2CYe,OnlyDiagonal) 
 TpYeCYeTpYeCYeml2 = Matmul2(Transpose(Ye),CYeTpYeCYeml2,OnlyDiagonal) 
 TpYumu2CYuTpYuCYu = Matmul2(Transpose(Yu),mu2CYuTpYuCYu,OnlyDiagonal) 
 TpYuCYumq2TpYuCYu = Matmul2(Transpose(Yu),CYumq2TpYuCYu,OnlyDiagonal) 
 TpYuCYuTpYumu2CYu = Matmul2(Transpose(Yu),CYuTpYumu2CYu,OnlyDiagonal) 
 TpYuCYuTpYuCYumq2 = Matmul2(Transpose(Yu),CYuTpYuCYumq2,OnlyDiagonal) 
 TpYw3mHw32CYw3TpYb3CYb3 = Matmul2(Transpose(Yw3),mHw32CYw3TpYb3CYb3,OnlyDiagonal) 
 TpYw3mHw32CYw3TpYw3CYw3 = Matmul2(Transpose(Yw3),mHw32CYw3TpYw3CYw3,OnlyDiagonal) 
 TpYw3CYw3ml2TpYb3CYb3 = Matmul2(Transpose(Yw3),CYw3ml2TpYb3CYb3,OnlyDiagonal) 
 TpYw3CYw3ml2TpYw3CYw3 = Matmul2(Transpose(Yw3),CYw3ml2TpYw3CYw3,OnlyDiagonal) 
 TpYw3CYw3TpYb3mHb32CYb3 = Matmul2(Transpose(Yw3),CYw3TpYb3mHb32CYb3,OnlyDiagonal) 
 TpYw3CYw3TpYb3CYb3ml2 = Matmul2(Transpose(Yw3),CYw3TpYb3CYb3ml2,OnlyDiagonal) 
 TpYw3CYw3TpYw3mHw32CYw3 = Matmul2(Transpose(Yw3),CYw3TpYw3mHw32CYw3,OnlyDiagonal) 
 TpYw3CYw3TpYw3CYw3ml2 = Matmul2(Transpose(Yw3),CYw3TpYw3CYw3ml2,OnlyDiagonal) 
 TpYx3mHxb32CYx3YdadjYd = Matmul2(Transpose(Yx3),mHxb32CYx3YdadjYd,OnlyDiagonal) 
 TpYx3mHxb32CYx3TpYx3CYx3 = Matmul2(Transpose(Yx3),mHxb32CYx3TpYx3CYx3,OnlyDiagonal) 
 TpYx3CYx3md2YdadjYd = Matmul2(Transpose(Yx3),CYx3md2YdadjYd,OnlyDiagonal) 
 TpYx3CYx3md2TpYx3CYx3 = Matmul2(Transpose(Yx3),CYx3md2TpYx3CYx3,OnlyDiagonal) 
 TpYx3CYx3TpYx3mHxb32CYx3 = Matmul2(Transpose(Yx3),CYx3TpYx3mHxb32CYx3,OnlyDiagonal) 
 TpYx3CYx3TpYx3CYx3md2 = Matmul2(Transpose(Yx3),CYx3TpYx3CYx3md2,OnlyDiagonal) 
 TpYx3CYx3TpYx3CYx3Yd = Matmul2(Transpose(Yx3),CYx3TpYx3CYx3Yd,OnlyDiagonal) 
 TpYx3CYx3TpYx3CYx3TYd = Matmul2(Transpose(Yx3),CYx3TpYx3CYx3TYd,OnlyDiagonal) 
 TpYx3CYx3TpTYx3CYx3Yd = Matmul2(Transpose(Yx3),CYx3TpTYx3CYx3Yd,OnlyDiagonal) 
 TpTYx3CYx3TpYx3CYx3Yd = Matmul2(Transpose(TYx3),CYx3TpYx3CYx3Yd,OnlyDiagonal) 
 TrmHg32 = cTrace(mHg32) 
 TrmHw32 = cTrace(mHw32) 
 TrTpYb3CTYb3 = cTrace(TpYb3CTYb3) 
 TrTpYdCTYd = cTrace(TpYdCTYd) 
 TrTpYeCTYe = cTrace(TpYeCTYe) 
 TrTpYuCTYu = cTrace(TpYuCTYu) 
 TrTpYw3CTYw3 = cTrace(TpYw3CTYw3) 
 TrTpYx3CTYx3 = cTrace(TpYx3CTYx3) 
 TrYb3adjYb3Yb3adjYb3 = cTrace(Yb3adjYb3Yb3adjYb3) 
 TrYb3adjYb3TYb3adjYb3 = cTrace(Yb3adjYb3TYb3adjYb3) 
 TrYb3adjYb3TYb3adjTYb3 = cTrace(Yb3adjYb3TYb3adjTYb3) 
 TrYb3adjYeYeadjYb3 = cTrace(Yb3adjYeYeadjYb3) 
 TrYb3adjYeTYeadjYb3 = cTrace(Yb3adjYeTYeadjYb3) 
 TrYb3adjYeTYeadjTYb3 = cTrace(Yb3adjYeTYeadjTYb3) 
 TrYb3adjYw3Yw3adjYb3 = cTrace(Yb3adjYw3Yw3adjYb3) 
 TrYb3adjYw3TYw3adjYb3 = cTrace(Yb3adjYw3TYw3adjYb3) 
 TrYb3adjYw3TYw3adjTYb3 = cTrace(Yb3adjYw3TYw3adjTYb3) 
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
 TrYdTpTYdCTYdadjYd = cTrace(YdTpTYdCTYdadjYd) 
 TrYdTpTYuCTYuadjYd = cTrace(YdTpTYuCTYuadjYd) 
 TrYeadjYb3Yb3adjYe = cTrace(YeadjYb3Yb3adjYe) 
 TrYeadjYb3TYb3adjYe = cTrace(YeadjYb3TYb3adjYe) 
 TrYeadjYb3TYb3adjTYe = cTrace(YeadjYb3TYb3adjTYe) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYeTYeadjYe = cTrace(YeadjYeTYeadjYe) 
 TrYeadjYeTYeadjTYe = cTrace(YeadjYeTYeadjTYe) 
 TrYeadjYw3Yw3adjYe = cTrace(YeadjYw3Yw3adjYe) 
 TrYeadjYw3TYw3adjYe = cTrace(YeadjYw3TYw3adjYe) 
 TrYeadjYw3TYw3adjTYe = cTrace(YeadjYw3TYw3adjTYe) 
 TrYeTpTYb3CTYb3adjYe = cTrace(YeTpTYb3CTYb3adjYe) 
 TrYeTpTYeCTYeadjYe = cTrace(YeTpTYeCTYeadjYe) 
 TrYeTpTYw3CTYw3adjYe = cTrace(YeTpTYw3CTYw3adjYe) 
 TrYuadjYdYdadjYu = cTrace(YuadjYdYdadjYu) 
 TrYuadjYdTYdadjYu = cTrace(YuadjYdTYdadjYu) 
 TrYuadjYdTYdadjTYu = cTrace(YuadjYdTYdadjTYu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYuadjYuTYuadjYu = cTrace(YuadjYuTYuadjYu) 
 TrYuadjYuTYuadjTYu = cTrace(YuadjYuTYuadjTYu) 
 TrYuTpTYdCTYdadjYu = cTrace(YuTpTYdCTYdadjYu) 
 TrYuTpTYuCTYuadjYu = cTrace(YuTpTYuCTYuadjYu) 
 TrYw3adjYb3Yb3adjYw3 = cTrace(Yw3adjYb3Yb3adjYw3) 
 TrYw3adjYb3TYb3adjYw3 = cTrace(Yw3adjYb3TYb3adjYw3) 
 TrYw3adjYb3TYb3adjTYw3 = cTrace(Yw3adjYb3TYb3adjTYw3) 
 TrYw3adjYeYeadjYw3 = cTrace(Yw3adjYeYeadjYw3) 
 TrYw3adjYeTYeadjYw3 = cTrace(Yw3adjYeTYeadjYw3) 
 TrYw3adjYeTYeadjTYw3 = cTrace(Yw3adjYeTYeadjTYw3) 
 TrYw3adjYw3Yw3adjYw3 = cTrace(Yw3adjYw3Yw3adjYw3) 
 TrYw3adjYw3TYw3adjYw3 = cTrace(Yw3adjYw3TYw3adjYw3) 
 TrYw3adjYw3TYw3adjTYw3 = cTrace(Yw3adjYw3TYw3adjTYw3) 
 TrYw3TpTYb3CTYb3adjYw3 = cTrace(Yw3TpTYb3CTYb3adjYw3) 
 TrYw3TpTYeCTYeadjYw3 = cTrace(Yw3TpTYeCTYeadjYw3) 
 TrYw3TpTYw3CTYw3adjYw3 = cTrace(Yw3TpTYw3CTYw3adjYw3) 
 TrYx3adjYx3Yx3adjYx3 = cTrace(Yx3adjYx3Yx3adjYx3) 
 TrYx3adjYx3TYx3adjYx3 = cTrace(Yx3adjYx3TYx3adjYx3) 
 TrYx3adjYx3TYx3adjTYx3 = cTrace(Yx3adjYx3TYx3adjTYx3) 
 TrYx3TpTYx3CTYx3adjYx3 = cTrace(Yx3TpTYx3CTYx3adjYx3) 
 TrCYdTpYdadjYx3TYx3 = cTrace(CYdTpYdadjYx3TYx3) 
 TrCTYdTpYdadjYx3TYx3 = cTrace(CTYdTpYdadjYx3TYx3) 
 TrTYb3adjYb3Yb3adjTYb3 = cTrace(TYb3adjYb3Yb3adjTYb3) 
 TrTYb3adjYeYeadjTYb3 = cTrace(TYb3adjYeYeadjTYb3) 
 TrTYb3adjYw3Yw3adjTYb3 = cTrace(TYb3adjYw3Yw3adjTYb3) 
 TrTYdadjYdYdadjTYd = cTrace(TYdadjYdYdadjTYd) 
 TrTYdadjYdTpYx3CYx3 = cTrace(TYdadjYdTpYx3CYx3) 
 TrTYdadjYdTpYx3CTYx3 = cTrace(TYdadjYdTpYx3CTYx3) 
 TrTYdadjYuYuadjTYd = cTrace(TYdadjYuYuadjTYd) 
 TrTYeadjYb3Yb3adjTYe = cTrace(TYeadjYb3Yb3adjTYe) 
 TrTYeadjYeYeadjTYe = cTrace(TYeadjYeYeadjTYe) 
 TrTYeadjYw3Yw3adjTYe = cTrace(TYeadjYw3Yw3adjTYe) 
 TrTYuadjYdYdadjTYu = cTrace(TYuadjYdYdadjTYu) 
 TrTYuadjYuYuadjTYu = cTrace(TYuadjYuYuadjTYu) 
 TrTYw3adjYb3Yb3adjTYw3 = cTrace(TYw3adjYb3Yb3adjTYw3) 
 TrTYw3adjYeYeadjTYw3 = cTrace(TYw3adjYeYeadjTYw3) 
 TrTYw3adjYw3Yw3adjTYw3 = cTrace(TYw3adjYw3Yw3adjTYw3) 
 TrTYx3adjYx3Yx3adjTYx3 = cTrace(TYx3adjYx3Yx3adjTYx3) 
 TrTpYdadjYdTpTYx3CTYx3 = cTrace(TpYdadjYdTpTYx3CTYx3) 
 TrTpYdadjYx3TYx3CTYd = cTrace(TpYdadjYx3TYx3CTYd) 
 TrTpYx3CYx3YdadjYd = cTrace(TpYx3CYx3YdadjYd) 
 TrTpYx3CYx3TYdadjYd = cTrace(TpYx3CYx3TYdadjYd) 
 TrTpTYdadjYdTpYx3CTYx3 = cTrace(TpTYdadjYdTpYx3CTYx3) 
 TrTpTYdadjYx3Yx3CTYd = cTrace(TpTYdadjYx3Yx3CTYd) 
 Trmd2YdadjYdTpYx3CYx3 = cTrace(md2YdadjYdTpYx3CYx3) 
 Trmd2TpYx3CYx3YdadjYd = cTrace(md2TpYx3CYx3YdadjYd) 
 Trmq2adjYdTpYx3CYx3Yd = cTrace(mq2adjYdTpYx3CYx3Yd) 
 TrYb3ml2adjYb3Yb3adjYb3 = cTrace(Yb3ml2adjYb3Yb3adjYb3) 
 TrYb3ml2adjYeYeadjYb3 = cTrace(Yb3ml2adjYeYeadjYb3) 
 TrYb3ml2adjYw3Yw3adjYb3 = cTrace(Yb3ml2adjYw3Yw3adjYb3) 
 TrYb3adjYb3mHb32Yb3adjYb3 = cTrace(Yb3adjYb3mHb32Yb3adjYb3) 
 TrYb3adjYb3Yb3ml2adjYb3 = cTrace(Yb3adjYb3Yb3ml2adjYb3) 
 TrYb3adjYeme2YeadjYb3 = cTrace(Yb3adjYeme2YeadjYb3) 
 TrYb3adjYeYeml2adjYb3 = cTrace(Yb3adjYeYeml2adjYb3) 
 TrYb3adjYw3mHw32Yw3adjYb3 = cTrace(Yb3adjYw3mHw32Yw3adjYb3) 
 TrYb3adjYw3Yw3ml2adjYb3 = cTrace(Yb3adjYw3Yw3ml2adjYb3) 
 TrYdmq2adjYdYdadjYd = cTrace(Ydmq2adjYdYdadjYd) 
 TrYdmq2adjYdTpYx3CYx3 = cTrace(Ydmq2adjYdTpYx3CYx3) 
 TrYdmq2adjYuYuadjYd = cTrace(Ydmq2adjYuYuadjYd) 
 TrYdadjYdmd2YdadjYd = cTrace(YdadjYdmd2YdadjYd) 
 TrYdadjYdmd2TpYx3CYx3 = cTrace(YdadjYdmd2TpYx3CYx3) 
 TrYdadjYdYdmq2adjYd = cTrace(YdadjYdYdmq2adjYd) 
 TrYdadjYumu2YuadjYd = cTrace(YdadjYumu2YuadjYd) 
 TrYdadjYuYumq2adjYd = cTrace(YdadjYuYumq2adjYd) 
 TrYeml2adjYb3Yb3adjYe = cTrace(Yeml2adjYb3Yb3adjYe) 
 TrYeml2adjYeYeadjYe = cTrace(Yeml2adjYeYeadjYe) 
 TrYeml2adjYw3Yw3adjYe = cTrace(Yeml2adjYw3Yw3adjYe) 
 TrYeadjYb3mHb32Yb3adjYe = cTrace(YeadjYb3mHb32Yb3adjYe) 
 TrYeadjYb3Yb3ml2adjYe = cTrace(YeadjYb3Yb3ml2adjYe) 
 TrYeadjYeme2YeadjYe = cTrace(YeadjYeme2YeadjYe) 
 TrYeadjYeYeml2adjYe = cTrace(YeadjYeYeml2adjYe) 
 TrYeadjYw3mHw32Yw3adjYe = cTrace(YeadjYw3mHw32Yw3adjYe) 
 TrYeadjYw3Yw3ml2adjYe = cTrace(YeadjYw3Yw3ml2adjYe) 
 TrYumq2adjYdYdadjYu = cTrace(Yumq2adjYdYdadjYu) 
 TrYumq2adjYuYuadjYu = cTrace(Yumq2adjYuYuadjYu) 
 TrYuadjYdmd2YdadjYu = cTrace(YuadjYdmd2YdadjYu) 
 TrYuadjYdYdmq2adjYu = cTrace(YuadjYdYdmq2adjYu) 
 TrYuadjYumu2YuadjYu = cTrace(YuadjYumu2YuadjYu) 
 TrYuadjYuYumq2adjYu = cTrace(YuadjYuYumq2adjYu) 
 TrYw3ml2adjYb3Yb3adjYw3 = cTrace(Yw3ml2adjYb3Yb3adjYw3) 
 TrYw3ml2adjYeYeadjYw3 = cTrace(Yw3ml2adjYeYeadjYw3) 
 TrYw3ml2adjYw3Yw3adjYw3 = cTrace(Yw3ml2adjYw3Yw3adjYw3) 
 TrYw3adjYb3mHb32Yb3adjYw3 = cTrace(Yw3adjYb3mHb32Yb3adjYw3) 
 TrYw3adjYb3Yb3ml2adjYw3 = cTrace(Yw3adjYb3Yb3ml2adjYw3) 
 TrYw3adjYeme2YeadjYw3 = cTrace(Yw3adjYeme2YeadjYw3) 
 TrYw3adjYeYeml2adjYw3 = cTrace(Yw3adjYeYeml2adjYw3) 
 TrYw3adjYw3mHw32Yw3adjYw3 = cTrace(Yw3adjYw3mHw32Yw3adjYw3) 
 TrYw3adjYw3Yw3ml2adjYw3 = cTrace(Yw3adjYw3Yw3ml2adjYw3) 
 TrYx3md2adjYx3Yx3adjYx3 = cTrace(Yx3md2adjYx3Yx3adjYx3) 
 TrYx3adjYx3mHxb32Yx3adjYx3 = cTrace(Yx3adjYx3mHxb32Yx3adjYx3) 
 TrYx3adjYx3Yx3md2adjYx3 = cTrace(Yx3adjYx3Yx3md2adjYx3) 
 TrCYdTpYdadjYx3mHxb32Yx3 = cTrace(CYdTpYdadjYx3mHxb32Yx3) 
 TrTpYx3mHxb32CYx3YdadjYd = cTrace(TpYx3mHxb32CYx3YdadjYd) 
 TrTpYx3CYx3md2YdadjYd = cTrace(TpYx3CYx3md2YdadjYd) 
 g1p4 =g1**4 
 g1p5 =g1**5 
 g2p4 =g2**4 
 g2p5 =g2**5 
 g3p4 =g3**4 
 g3p5 =g3**5 
End If 
 
 
Tr1(1) = -3._dp*mHd2/5._dp + 3._dp*mHu2/5._dp + 3._dp*Trmd2/5._dp + 3._dp*Trme2/5._dp +& 
&  3._dp*TrmHx32 - 3._dp*TrmHxb32 - 3._dp*Trml2/5._dp + 3._dp*Trmq2/5._dp -              & 
&  6._dp*Trmu2/5._dp

If (TwoLoopRGE) Then 
Tr2(1) = 3._dp*mHd2/10._dp + 3._dp*mHu2/10._dp + Trmd2/5._dp + 3._dp*Trme2/5._dp +    & 
&  5._dp*TrmHx32/2._dp + 5._dp*TrmHxb32/2._dp + 3._dp*Trml2/10._dp + Trmq2/10._dp +      & 
&  4._dp*Trmu2/5._dp

Tr3(1) = -9._dp*g1p2*mHd2/100._dp - 9._dp*g2p2*mHd2/20._dp + 9._dp*g1p2*mHu2/100._dp +& 
&  9._dp*g2p2*mHu2/20._dp + (g1p2*Trmd2)/25._dp + 4._dp*g3p2*Trmd2/5._dp -               & 
&  3._dp*Trmd2YdadjYd/5._dp + 9._dp*g1p2*Trme2/25._dp - 3._dp*Trme2YeadjYe/5._dp +       & 
&  5._dp*g1p2*TrmHx32/4._dp + 9._dp*g2p2*TrmHx32/4._dp + 4._dp*g3p2*TrmHx32 -            & 
&  5._dp*g1p2*TrmHxb32/4._dp - 9._dp*g2p2*TrmHxb32/4._dp - 4._dp*g3p2*TrmHxb32 +         & 
&  3._dp*TrmHxb32Yx3adjYx3/2._dp - 9._dp*g1p2*Trml2/100._dp - 9._dp*g2p2*Trml2/20._dp +  & 
&  (g1p2*Trmq2)/100._dp + 9._dp*g2p2*Trmq2/20._dp + 4._dp*g3p2*Trmq2/5._dp -             & 
&  8._dp*g1p2*Trmu2/25._dp - 8._dp*g3p2*Trmu2/5._dp + 6._dp*Trmu2YuadjYu/5._dp -         & 
&  9._dp*mHu2*TrYb3adjYb3/100._dp + 9._dp*TrYb3ml2adjYb3/100._dp + 9._dp*mHd2*TrYdadjYd/10._dp -& 
&  3._dp*TrYdmq2adjYd/10._dp + 3._dp*mHd2*TrYeadjYe/10._dp + 3._dp*TrYeml2adjYe/10._dp - & 
&  9._dp*mHu2*TrYuadjYu/10._dp - 3._dp*TrYumq2adjYu/10._dp - 9._dp*mHu2*TrYw3adjYw3/20._dp +& 
&  9._dp*TrYw3ml2adjYw3/20._dp - 9._dp*mHu2*TrYx3adjYx3/10._dp - 3._dp*TrYx3md2adjYx3/5._dp

Tr2(2) = mHd2/2._dp + mHu2/2._dp + 3._dp*TrmHw32 + 3._dp*TrmHx32/2._dp +              & 
&  3._dp*TrmHxb32/2._dp + Trml2/2._dp + 3._dp*Trmq2/2._dp

Tr2(3) = Trmd2/2._dp + 8._dp*TrmHg32 + TrmHx32 + TrmHxb32 + Trmq2 + Trmu2/2._dp

End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11  = 33._dp*g1p3/5._dp + 5._dp*g1p3*NGHx3/2._dp + 5._dp*g1p3*NGHxb3/2._dp

 
 
If (TwoLoopRGE) Then 
betag12 = 199._dp*g1p5/25._dp + 27._dp*g1p3*g2p2/5._dp + 88._dp*g1p3*g3p2/5._dp -               & 
&  9._dp*g1p3*TrYb3adjYb3/25._dp - 14._dp*g1p3*TrYdadjYd/5._dp - 18._dp*g1p3*TrYeadjYe/5._dp -& 
&  26._dp*g1p3*TrYuadjYu/5._dp - 12._dp*g1p3*TrYw3adjYw3/5._dp - 38._dp*g1p3*TrYx3adjYx3/5._dp +& 
&  25._dp*g1p5*NGHx3/6._dp + 15._dp*g1p3*g2p2*NGHx3/2._dp +& 
&  40._dp*g1p3*g3p2*NGHx3/3._dp + 25._dp*g1p5*NGHxb3/6._dp +& 
&  15._dp*g1p3*g2p2*NGHxb3/2._dp + 40._dp*g1p3*g3p2*NGHxb3/3._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21  = g2p3 + 2._dp*g2p3*NGHw3 + 3._dp*g2p3*NGHx3/2._dp +& 
&  3._dp*g2p3*NGHxb3/2._dp

 
 
If (TwoLoopRGE) Then 
betag22 = 9._dp*g1p2*g2p3/5._dp + 25._dp*g2p5 + 24._dp*g2p3*g3p2 - 3._dp*g2p3*TrYb3adjYb3/5._dp -& 
&  6._dp*g2p3*TrYdadjYd - 2._dp*g2p3*TrYeadjYe - 6._dp*g2p3*TrYuadjYu - 28._dp*g2p3*TrYw3adjYw3/3._dp -& 
&  6._dp*g2p3*TrYx3adjYx3 + 24._dp*g2p5*NGHw3 + 5._dp*g1p2*g2p3*NGHx3/2._dp +& 
&  21._dp*g2p5*NGHx3/2._dp + 8._dp*g2p3*g3p2*NGHx3 +       & 
&  5._dp*g1p2*g2p3*NGHxb3/2._dp + 21._dp*g2p5*NGHxb3/2._dp +& 
&  8._dp*g2p3*g3p2*NGHxb3

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31  = -3._dp*g3p3 + 3._dp*g3p3*NGHg3 + g3p3*NGHx3 +& 
&  g3p3*NGHxb3

 
 
If (TwoLoopRGE) Then 
betag32 = 11._dp*g1p2*g3p3/5._dp + 9._dp*g2p2*g3p3 + 14._dp*g3p5 - 4._dp*g3p3*TrYdadjYd -       & 
&  4._dp*g3p3*TrYuadjYu - 4._dp*g3p3*TrYx3adjYx3 + 54._dp*g3p5*NGHg3 +    & 
&  5._dp*g1p2*g3p3*NGHx3/3._dp + 3._dp*g2p2*g3p3*NGHx3 +   & 
&  34._dp*g3p5*NGHx3/3._dp + 5._dp*g1p2*g3p3*NGHxb3/3._dp +& 
&  3._dp*g2p2*g3p3*NGHxb3 + 34._dp*g3p5*NGHxb3/3._dp

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
betaYu1  = -13._dp*g1p2*Yu/15._dp - 3._dp*g2p2*Yu - 16._dp*g3p2*Yu/3._dp +            & 
&  3._dp*TrYb3adjYb3*Yu/10._dp + 3._dp*TrYuadjYu*Yu + 3._dp*TrYw3adjYw3*Yu/2._dp +       & 
&  3._dp*TrYx3adjYx3*Yu + YuadjYdYd + 3._dp*YuadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYu2 = 2743._dp*g1p4*Yu/450._dp + g1p2*g2p2*Yu + 15._dp*g2p4*Yu/2._dp + 136._dp*g1p2*g3p2*Yu/45._dp +& 
&  8._dp*g2p2*g3p2*Yu - 16._dp*g3p4*Yu/9._dp - 27._dp*TrYb3adjYb3Yb3adjYb3*Yu/100._dp -  & 
&  3._dp*TrYb3adjYeYeadjYb3*Yu/10._dp - 3._dp*TrYb3adjYw3Yw3adjYb3*Yu/4._dp -            & 
&  6._dp*TrYdadjYdTpYx3CYx3*Yu - 3._dp*TrYuadjYdYdadjYu*Yu + 4._dp*g1p2*TrYuadjYu*Yu/5._dp +& 
&  16._dp*g3p2*TrYuadjYu*Yu - 9._dp*TrYuadjYuYuadjYu*Yu - 27._dp*TrYw3adjYb3Yb3adjYw3*Yu/40._dp -& 
&  3._dp*TrYw3adjYeYeadjYw3*Yu/2._dp + 6._dp*g2p2*TrYw3adjYw3*Yu - 15._dp*TrYw3adjYw3Yw3adjYw3*Yu/4._dp +& 
&  2._dp*g1p2*TrYx3adjYx3*Yu + 16._dp*g3p2*TrYx3adjYx3*Yu - 9._dp*TrYx3adjYx3Yx3adjYx3*Yu -& 
&  2._dp*YuadjYdTpYx3CYx3Yd + 2._dp*g1p2*YuadjYdYd/5._dp - 3._dp*TrYdadjYd*YuadjYdYd -   & 
&  TrYeadjYe*YuadjYdYd - 2._dp*YuadjYdYdadjYdYd - 2._dp*YuadjYdYdadjYuYu +               & 
&  2._dp*g1p2*YuadjYuYu/5._dp + 6._dp*g2p2*YuadjYuYu - 9._dp*TrYb3adjYb3*YuadjYuYu/10._dp -& 
&  9._dp*TrYuadjYu*YuadjYuYu - 9._dp*TrYw3adjYw3*YuadjYuYu/2._dp - 9._dp*TrYx3adjYx3*YuadjYuYu  
betaYu2 =  betaYu2- 4._dp*YuadjYuYuadjYuYu + 16._dp*g3p4*Yu*NGHg3 + 6._dp*g2p4*Yu*NGHw3 +& 
&  13._dp*g1p4*Yu*NGHx3/6._dp + 9._dp*g2p4*Yu*NGHx3/2._dp +& 
&  16._dp*g3p4*Yu*NGHx3/3._dp + 13._dp*g1p4*Yu*NGHxb3/6._dp +& 
&  9._dp*g2p4*Yu*NGHxb3/2._dp + 16._dp*g3p4*Yu*NGHxb3/3._dp

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
betaYd1  = 2._dp*TpYx3CYx3Yd - 7._dp*g1p2*Yd/15._dp - 3._dp*g2p2*Yd - 16._dp*g3p2*Yd/3._dp +& 
&  3._dp*TrYdadjYd*Yd + TrYeadjYe*Yd + 3._dp*YdadjYdYd + YdadjYuYu

 
 
If (TwoLoopRGE) Then 
betaYd2 = -2._dp*TpYx3CYx3TpYx3CYx3Yd + 2._dp*g1p2*TpYx3CYx3Yd + 6._dp*g2p2*TpYx3CYx3Yd -       & 
&  3._dp*TpYx3CYx3Yd*TrYb3adjYb3/5._dp - 6._dp*TpYx3CYx3Yd*TrYuadjYu - 3._dp*TpYx3CYx3Yd*TrYw3adjYw3 -& 
&  6._dp*TpYx3CYx3Yd*TrYx3adjYx3 + 287._dp*g1p4*Yd/90._dp + g1p2*g2p2*Yd +               & 
&  15._dp*g2p4*Yd/2._dp + 8._dp*g1p2*g3p2*Yd/9._dp + 8._dp*g2p2*g3p2*Yd - 16._dp*g3p4*Yd/9._dp -& 
&  6._dp*TrTpYx3CYx3YdadjYd*Yd - 2._dp*g1p2*TrYdadjYd*Yd/5._dp + 16._dp*g3p2*TrYdadjYd*Yd -& 
&  9._dp*TrYdadjYdYdadjYd*Yd - 3._dp*TrYdadjYuYuadjYd*Yd - 3._dp*TrYeadjYb3Yb3adjYe*Yd/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*Yd/5._dp - 3._dp*TrYeadjYeYeadjYe*Yd - 3._dp*TrYeadjYw3Yw3adjYe*Yd/2._dp -& 
&  2._dp*YdadjYdTpYx3CYx3Yd + 4._dp*g1p2*YdadjYdYd/5._dp + 6._dp*g2p2*YdadjYdYd -        & 
&  9._dp*TrYdadjYd*YdadjYdYd - 3._dp*TrYeadjYe*YdadjYdYd - 4._dp*YdadjYdYdadjYdYd +      & 
&  4._dp*g1p2*YdadjYuYu/5._dp - 3._dp*TrYb3adjYb3*YdadjYuYu/10._dp - 3._dp*TrYuadjYu*YdadjYuYu -& 
&  3._dp*TrYw3adjYw3*YdadjYuYu/2._dp - 3._dp*TrYx3adjYx3*YdadjYuYu - 2._dp*YdadjYuYuadjYdYd  
betaYd2 =  betaYd2- 2._dp*YdadjYuYuadjYuYu + 16._dp*g3p4*Yd*NGHg3 + 6._dp*g2p4*Yd*NGHw3 +& 
&  7._dp*g1p4*Yd*NGHx3/6._dp + 9._dp*g2p4*Yd*NGHx3/2._dp + & 
&  16._dp*g3p4*Yd*NGHx3/3._dp + 7._dp*g1p4*Yd*NGHxb3/6._dp +& 
&  9._dp*g2p4*Yd*NGHxb3/2._dp + 16._dp*g3p4*Yd*NGHxb3/3._dp

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
betaYe1  = -9._dp*g1p2*Ye/5._dp - 3._dp*g2p2*Ye + 3._dp*TrYdadjYd*Ye + TrYeadjYe*Ye + & 
&  3._dp*YeadjYb3Yb3/10._dp + 3._dp*YeadjYeYe + 3._dp*YeadjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYe2 = 27._dp*g1p4*Ye/2._dp + 9._dp*g1p2*g2p2*Ye/5._dp + 15._dp*g2p4*Ye/2._dp -              & 
&  6._dp*TrTpYx3CYx3YdadjYd*Ye - 2._dp*g1p2*TrYdadjYd*Ye/5._dp + 16._dp*g3p2*TrYdadjYd*Ye -& 
&  9._dp*TrYdadjYdYdadjYd*Ye - 3._dp*TrYdadjYuYuadjYd*Ye - 3._dp*TrYeadjYb3Yb3adjYe*Ye/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*Ye/5._dp - 3._dp*TrYeadjYeYeadjYe*Ye - 3._dp*TrYeadjYw3Yw3adjYe*Ye/2._dp -& 
&  9._dp*TrYb3adjYb3*YeadjYb3Yb3/100._dp - 9._dp*TrYuadjYu*YeadjYb3Yb3/10._dp -          & 
&  9._dp*TrYw3adjYw3*YeadjYb3Yb3/20._dp - 9._dp*TrYx3adjYx3*YeadjYb3Yb3/10._dp -         & 
&  9._dp*YeadjYb3Yb3adjYb3Yb3/50._dp - 3._dp*YeadjYb3Yb3adjYeYe/5._dp - 3._dp*YeadjYb3Yb3adjYw3Yw3/10._dp +& 
&  6._dp*g2p2*YeadjYeYe - 9._dp*TrYdadjYd*YeadjYeYe - 3._dp*TrYeadjYe*YeadjYeYe -        & 
&  4._dp*YeadjYeYeadjYeYe + 6._dp*g2p2*YeadjYw3Yw3 - 9._dp*TrYb3adjYb3*YeadjYw3Yw3/20._dp -& 
&  9._dp*TrYuadjYu*YeadjYw3Yw3/2._dp - 9._dp*TrYw3adjYw3*YeadjYw3Yw3/4._dp -             & 
&  9._dp*TrYx3adjYx3*YeadjYw3Yw3/2._dp - 9._dp*YeadjYw3Yw3adjYb3Yb3/40._dp  
betaYe2 =  betaYe2- 3._dp*YeadjYw3Yw3adjYeYe - 3._dp*YeadjYw3Yw3adjYw3Yw3/2._dp + 6._dp*g2p4*Ye*NGHw3 +& 
&  9._dp*g1p4*Ye*NGHx3/2._dp + 9._dp*g2p4*Ye*NGHx3/2._dp + & 
&  9._dp*g1p4*Ye*NGHxb3/2._dp + 9._dp*g2p4*Ye*NGHxb3/2._dp

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yb3 
!-------------------- 
 
betaYb31  = -3._dp*g1p2*Yb3/5._dp - 3._dp*g2p2*Yb3 + 3._dp*TrYb3adjYb3*Yb3/10._dp +   & 
&  3._dp*TrYuadjYu*Yb3 + 3._dp*TrYw3adjYw3*Yb3/2._dp + 3._dp*TrYx3adjYx3*Yb3 +           & 
&  9._dp*Yb3adjYb3Yb3/10._dp + Yb3adjYeYe + 3._dp*Yb3adjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYb32 = 207._dp*g1p4*Yb3/50._dp + 9._dp*g1p2*g2p2*Yb3/5._dp + 15._dp*g2p4*Yb3/2._dp -         & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yb3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yb3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yb3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yb3 - 3._dp*TrYuadjYdYdadjYu*Yb3 +& 
&  4._dp*g1p2*TrYuadjYu*Yb3/5._dp + 16._dp*g3p2*TrYuadjYu*Yb3 - 9._dp*TrYuadjYuYuadjYu*Yb3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yb3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yb3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yb3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yb3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yb3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yb3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yb3 + 9._dp*g1p2*Yb3adjYb3Yb3/25._dp +& 
&  9._dp*g2p2*Yb3adjYb3Yb3/5._dp - 27._dp*TrYb3adjYb3*Yb3adjYb3Yb3/100._dp -             & 
&  27._dp*TrYuadjYu*Yb3adjYb3Yb3/10._dp - 27._dp*TrYw3adjYw3*Yb3adjYb3Yb3/20._dp -       & 
&  27._dp*TrYx3adjYx3*Yb3adjYb3Yb3/10._dp - 9._dp*Yb3adjYb3Yb3adjYb3Yb3/25._dp -         & 
&  3._dp*Yb3adjYb3Yb3adjYw3Yw3/10._dp + 6._dp*g1p2*Yb3adjYeYe/5._dp - 3._dp*TrYdadjYd*Yb3adjYeYe  
betaYb32 =  betaYb32- TrYeadjYe*Yb3adjYeYe - 3._dp*Yb3adjYeYeadjYb3Yb3/5._dp - 2._dp*Yb3adjYeYeadjYeYe +    & 
&  6._dp*g2p2*Yb3adjYw3Yw3 - 9._dp*TrYb3adjYb3*Yb3adjYw3Yw3/20._dp - 9._dp*TrYuadjYu*Yb3adjYw3Yw3/2._dp -& 
&  9._dp*TrYw3adjYw3*Yb3adjYw3Yw3/4._dp - 9._dp*TrYx3adjYx3*Yb3adjYw3Yw3/2._dp -         & 
&  9._dp*Yb3adjYw3Yw3adjYb3Yb3/8._dp - 3._dp*Yb3adjYw3Yw3adjYw3Yw3/2._dp +               & 
&  6._dp*g2p4*Yb3*NGHw3 + 3._dp*g1p4*Yb3*NGHx3/2._dp +     & 
&  9._dp*g2p4*Yb3*NGHx3/2._dp + 3._dp*g1p4*Yb3*NGHxb3/2._dp +& 
&  9._dp*g2p4*Yb3*NGHxb3/2._dp

 
DYb3 = oo16pi2*( betaYb31 + oo16pi2 * betaYb32 ) 

 
Else 
DYb3 = oo16pi2* betaYb31 
End If 
 
 
!-------------------- 
! Yw3 
!-------------------- 
 
betaYw31  = -3._dp*g1p2*Yw3/5._dp - 7._dp*g2p2*Yw3 + 3._dp*TrYb3adjYb3*Yw3/10._dp +   & 
&  3._dp*TrYuadjYu*Yw3 + 3._dp*TrYw3adjYw3*Yw3/2._dp + 3._dp*TrYx3adjYx3*Yw3 +           & 
&  3._dp*Yw3adjYb3Yb3/10._dp + Yw3adjYeYe + 5._dp*Yw3adjYw3Yw3/2._dp

 
 
If (TwoLoopRGE) Then 
betaYw32 = 207._dp*g1p4*Yw3/50._dp + 9._dp*g1p2*g2p2*Yw3/5._dp + 55._dp*g2p4*Yw3/2._dp -         & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yw3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yw3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yw3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yw3 - 3._dp*TrYuadjYdYdadjYu*Yw3 +& 
&  4._dp*g1p2*TrYuadjYu*Yw3/5._dp + 16._dp*g3p2*TrYuadjYu*Yw3 - 9._dp*TrYuadjYuYuadjYu*Yw3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yw3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yw3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yw3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yw3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yw3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yw3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yw3 - 9._dp*TrYb3adjYb3*Yw3adjYb3Yb3/100._dp -& 
&  9._dp*TrYuadjYu*Yw3adjYb3Yb3/10._dp - 9._dp*TrYw3adjYw3*Yw3adjYb3Yb3/20._dp -         & 
&  9._dp*TrYx3adjYx3*Yw3adjYb3Yb3/10._dp - 9._dp*Yw3adjYb3Yb3adjYb3Yb3/50._dp -          & 
&  3._dp*Yw3adjYb3Yb3adjYw3Yw3/5._dp + 6._dp*g1p2*Yw3adjYeYe/5._dp - 3._dp*TrYdadjYd*Yw3adjYeYe -& 
&  TrYeadjYe*Yw3adjYeYe - 2._dp*Yw3adjYeYeadjYeYe - Yw3adjYeYeadjYw3Yw3 + 3._dp*g1p2*Yw3adjYw3Yw3/5._dp  
betaYw32 =  betaYw32+ 5._dp*g2p2*Yw3adjYw3Yw3 - 3._dp*TrYb3adjYb3*Yw3adjYw3Yw3/4._dp - 15._dp*TrYuadjYu*Yw3adjYw3Yw3/2._dp -& 
&  15._dp*TrYw3adjYw3*Yw3adjYw3Yw3/4._dp - 15._dp*TrYx3adjYx3*Yw3adjYw3Yw3/2._dp -       & 
&  9._dp*Yw3adjYw3Yw3adjYb3Yb3/40._dp - 3._dp*Yw3adjYw3Yw3adjYw3Yw3 + 14._dp*g2p4*Yw3*NGHw3 +& 
&  3._dp*g1p4*Yw3*NGHx3/2._dp + 21._dp*g2p4*Yw3*NGHx3/2._dp +& 
&  3._dp*g1p4*Yw3*NGHxb3/2._dp + 21._dp*g2p4*Yw3*NGHxb3/2._dp

 
DYw3 = oo16pi2*( betaYw31 + oo16pi2 * betaYw32 ) 

 
Else 
DYw3 = oo16pi2* betaYw31 
End If 
 
 
!-------------------- 
! Yx3 
!-------------------- 
 
betaYx31  = -19._dp*g1p2*Yx3/15._dp - 3._dp*g2p2*Yx3 - 16._dp*g3p2*Yx3/3._dp +        & 
&  3._dp*TrYb3adjYb3*Yx3/10._dp + 3._dp*TrYuadjYu*Yx3 + 3._dp*TrYw3adjYw3*Yx3/2._dp +    & 
&  3._dp*TrYx3adjYx3*Yx3 + 3._dp*Yx3adjYx3Yx3 + 2._dp*Yx3CYdTpYd

 
 
If (TwoLoopRGE) Then 
betaYx32 = 4123._dp*g1p4*Yx3/450._dp + 17._dp*g1p2*g2p2*Yx3/5._dp + 15._dp*g2p4*Yx3/2._dp +      & 
&  232._dp*g1p2*g3p2*Yx3/45._dp + 8._dp*g2p2*g3p2*Yx3 - 16._dp*g3p4*Yx3/9._dp -          & 
&  27._dp*TrYb3adjYb3Yb3adjYb3*Yx3/100._dp - 3._dp*TrYb3adjYeYeadjYb3*Yx3/10._dp -       & 
&  3._dp*TrYb3adjYw3Yw3adjYb3*Yx3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3*Yx3 - 3._dp*TrYuadjYdYdadjYu*Yx3 +& 
&  4._dp*g1p2*TrYuadjYu*Yx3/5._dp + 16._dp*g3p2*TrYuadjYu*Yx3 - 9._dp*TrYuadjYuYuadjYu*Yx3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*Yx3/40._dp - 3._dp*TrYw3adjYeYeadjYw3*Yx3/2._dp +         & 
&  6._dp*g2p2*TrYw3adjYw3*Yx3 - 15._dp*TrYw3adjYw3Yw3adjYw3*Yx3/4._dp + 2._dp*g1p2*TrYx3adjYx3*Yx3 +& 
&  16._dp*g3p2*TrYx3adjYx3*Yx3 - 9._dp*TrYx3adjYx3Yx3adjYx3*Yx3 + 8._dp*g1p2*Yx3adjYx3Yx3/5._dp +& 
&  6._dp*g2p2*Yx3adjYx3Yx3 - 9._dp*TrYb3adjYb3*Yx3adjYx3Yx3/10._dp - 9._dp*TrYuadjYu*Yx3adjYx3Yx3 -& 
&  9._dp*TrYw3adjYw3*Yx3adjYx3Yx3/2._dp - 9._dp*TrYx3adjYx3*Yx3adjYx3Yx3 -               & 
&  4._dp*Yx3adjYx3Yx3adjYx3Yx3 + 2._dp*g1p2*Yx3CYdTpYd/5._dp + 6._dp*g2p2*Yx3CYdTpYd  
betaYx32 =  betaYx32- 6._dp*TrYdadjYd*Yx3CYdTpYd - 2._dp*TrYeadjYe*Yx3CYdTpYd - 2._dp*Yx3CYdTpYdadjYx3Yx3 - & 
&  2._dp*Yx3CYdTpYdCYdTpYd - 2._dp*Yx3CYdTpYuCYuTpYd + 16._dp*g3p4*Yx3*NGHg3 +& 
&  6._dp*g2p4*Yx3*NGHw3 + 19._dp*g1p4*Yx3*NGHx3/6._dp +    & 
&  9._dp*g2p4*Yx3*NGHx3/2._dp + 16._dp*g3p4*Yx3*NGHx3/3._dp +& 
&  19._dp*g1p4*Yx3*NGHxb3/6._dp + 9._dp*g2p4*Yx3*NGHxb3/2._dp +& 
&  16._dp*g3p4*Yx3*NGHxb3/3._dp

 
DYx3 = oo16pi2*( betaYx31 + oo16pi2 * betaYx32 ) 

 
Else 
DYx3 = oo16pi2* betaYx31 
End If 
 
 
!-------------------- 
! mue 
!-------------------- 
 
betamue1  = -3._dp*g1p2*mue/5._dp - 3._dp*g2p2*mue + 3._dp*TrYb3adjYb3*mue/10._dp +   & 
&  3._dp*TrYdadjYd*mue + TrYeadjYe*mue + 3._dp*TrYuadjYu*mue + 3._dp*TrYw3adjYw3*mue/2._dp +& 
&  3._dp*TrYx3adjYx3*mue

 
 
If (TwoLoopRGE) Then 
betamue2 = 207._dp*g1p4*mue/50._dp + 9._dp*g1p2*g2p2*mue/5._dp + 15._dp*g2p4*mue/2._dp -         & 
&  6._dp*TrTpYx3CYx3YdadjYd*mue - 27._dp*TrYb3adjYb3Yb3adjYb3*mue/100._dp -              & 
&  3._dp*TrYb3adjYeYeadjYb3*mue/10._dp - 3._dp*TrYb3adjYw3Yw3adjYb3*mue/4._dp -          & 
&  2._dp*g1p2*TrYdadjYd*mue/5._dp + 16._dp*g3p2*TrYdadjYd*mue - 6._dp*TrYdadjYdTpYx3CYx3*mue -& 
&  9._dp*TrYdadjYdYdadjYd*mue - 3._dp*TrYdadjYuYuadjYd*mue - 3._dp*TrYeadjYb3Yb3adjYe*mue/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*mue/5._dp - 3._dp*TrYeadjYeYeadjYe*mue - 3._dp*TrYeadjYw3Yw3adjYe*mue/2._dp -& 
&  3._dp*TrYuadjYdYdadjYu*mue + 4._dp*g1p2*TrYuadjYu*mue/5._dp + 16._dp*g3p2*TrYuadjYu*mue -& 
&  9._dp*TrYuadjYuYuadjYu*mue - 27._dp*TrYw3adjYb3Yb3adjYw3*mue/40._dp - 3._dp*TrYw3adjYeYeadjYw3*mue/2._dp +& 
&  6._dp*g2p2*TrYw3adjYw3*mue - 15._dp*TrYw3adjYw3Yw3adjYw3*mue/4._dp + 2._dp*g1p2*TrYx3adjYx3*mue +& 
&  16._dp*g3p2*TrYx3adjYx3*mue - 9._dp*TrYx3adjYx3Yx3adjYx3*mue + 6._dp*g2p4*mue*NGHw3 +& 
&  3._dp*g1p4*mue*NGHx3/2._dp + 9._dp*g2p4*mue*NGHx3/2._dp  
betamue2 =  betamue2+ 3._dp*g1p4*mue*NGHxb3/2._dp + 9._dp*g2p4*mue*NGHxb3/2._dp

 
Dmue = oo16pi2*( betamue1 + oo16pi2 * betamue2 ) 

 
Else 
Dmue = oo16pi2* betamue1 
End If 
 
 
!-------------------- 
! MXM3 
!-------------------- 
 
betaMXM31  = -5._dp*g1p2*MXM3/3._dp - 3._dp*g2p2*MXM3 - 16._dp*g3p2*MXM3/3._dp +      & 
&  MXM3CYx3TpYx3

 
 
If (TwoLoopRGE) Then 
betaMXM32 = 223._dp*g1p4*MXM3/18._dp + 5._dp*g1p2*g2p2*MXM3 + 15._dp*g2p4*MXM3/2._dp +            & 
&  80._dp*g1p2*g3p2*MXM3/9._dp + 16._dp*g2p2*g3p2*MXM3 - 16._dp*g3p4*MXM3/9._dp -        & 
&  2._dp*g1p2*MXM3CYx3TpYx3/5._dp - 2._dp*MXM3CYx3TpYx3CYx3TpYx3 - 2._dp*MXM3CYx3YdadjYdTpYx3 -& 
&  3._dp*MXM3CYx3TpYx3*TrYb3adjYb3/10._dp - 3._dp*MXM3CYx3TpYx3*TrYuadjYu -              & 
&  3._dp*MXM3CYx3TpYx3*TrYw3adjYw3/2._dp - 3._dp*MXM3CYx3TpYx3*TrYx3adjYx3 +             & 
&  16._dp*g3p4*MXM3*NGHg3 + 6._dp*g2p4*MXM3*NGHw3 +        & 
&  25._dp*g1p4*MXM3*NGHx3/6._dp + 9._dp*g2p4*MXM3*NGHx3/2._dp +& 
&  16._dp*g3p4*MXM3*NGHx3/3._dp + 25._dp*g1p4*MXM3*NGHxb3/6._dp +& 
&  9._dp*g2p4*MXM3*NGHxb3/2._dp + 16._dp*g3p4*MXM3*NGHxb3/3._dp

 
DMXM3 = oo16pi2*( betaMXM31 + oo16pi2 * betaMXM32 ) 

 
Else 
DMXM3 = oo16pi2* betaMXM31 
End If 
 
 
!-------------------- 
! MWM3 
!-------------------- 
 
betaMWM31  = -8._dp*g2p2*MWM3 + MWM3CYw3TpYw3 + Yw3adjYw3MWM3

 
 
If (TwoLoopRGE) Then 
betaMWM32 = 40._dp*g2p4*MWM3 - 3._dp*MWM3CYw3TpYb3CYb3TpYw3/10._dp - MWM3CYw3TpYeCYeTpYw3 +       & 
&  3._dp*g1p2*MWM3CYw3TpYw3/5._dp - g2p2*MWM3CYw3TpYw3 - 3._dp*MWM3CYw3TpYw3CYw3TpYw3/2._dp -& 
&  3._dp*MWM3CYw3TpYw3*TrYb3adjYb3/10._dp - 3._dp*MWM3CYw3TpYw3*TrYuadjYu -              & 
&  3._dp*MWM3CYw3TpYw3*TrYw3adjYw3/2._dp - 3._dp*MWM3CYw3TpYw3*TrYx3adjYx3 -             & 
&  3._dp*Yw3adjYb3Yb3adjYw3MWM3/10._dp - Yw3adjYeYeadjYw3MWM3 + 3._dp*g1p2*Yw3adjYw3MWM3/5._dp -& 
&  g2p2*Yw3adjYw3MWM3 - 3._dp*TrYb3adjYb3*Yw3adjYw3MWM3/10._dp - 3._dp*TrYuadjYu*Yw3adjYw3MWM3 -& 
&  3._dp*TrYw3adjYw3*Yw3adjYw3MWM3/2._dp - 3._dp*TrYx3adjYx3*Yw3adjYw3MWM3 -             & 
&  3._dp*Yw3adjYw3Yw3adjYw3MWM3/2._dp + 16._dp*g2p4*MWM3*NGHw3 +          & 
&  12._dp*g2p4*MWM3*NGHx3 + 12._dp*g2p4*MWM3*NGHxb3

 
DMWM3 = oo16pi2*( betaMWM31 + oo16pi2 * betaMWM32 ) 

 
Else 
DMWM3 = oo16pi2* betaMWM31 
End If 
 
 
!-------------------- 
! MGM3 
!-------------------- 
 
betaMGM31  = -12._dp*g3p2*MGM3

 
 
If (TwoLoopRGE) Then 
betaMGM32 = 36._dp*g3p4*MGM3 + 36._dp*g3p4*MGM3*NGHg3 + 12._dp*g3p4*MGM3*NGHx3 +& 
&  12._dp*g3p4*MGM3*NGHxb3

 
DMGM3 = oo16pi2*( betaMGM31 + oo16pi2 * betaMGM32 ) 

 
Else 
DMGM3 = oo16pi2* betaMGM31 
End If 
 
 
!-------------------- 
! MBM3 
!-------------------- 
 
betaMBM31  = 3._dp*MBM3CYb3TpYb3/5._dp + 3._dp*Yb3adjYb3MBM3/5._dp

 
 
If (TwoLoopRGE) Then 
betaMBM32 = 9._dp*g1p2*MBM3CYb3TpYb3/25._dp + 9._dp*g2p2*MBM3CYb3TpYb3/5._dp - 9._dp*MBM3CYb3TpYb3CYb3TpYb3/50._dp -& 
&  3._dp*MBM3CYb3TpYeCYeTpYb3/5._dp - 9._dp*MBM3CYb3TpYw3CYw3TpYb3/10._dp -              & 
&  9._dp*MBM3CYb3TpYb3*TrYb3adjYb3/50._dp - 9._dp*MBM3CYb3TpYb3*TrYuadjYu/5._dp -        & 
&  9._dp*MBM3CYb3TpYb3*TrYw3adjYw3/10._dp - 9._dp*MBM3CYb3TpYb3*TrYx3adjYx3/5._dp +      & 
&  9._dp*g1p2*Yb3adjYb3MBM3/25._dp + 9._dp*g2p2*Yb3adjYb3MBM3/5._dp - 9._dp*TrYb3adjYb3*Yb3adjYb3MBM3/50._dp -& 
&  9._dp*TrYuadjYu*Yb3adjYb3MBM3/5._dp - 9._dp*TrYw3adjYw3*Yb3adjYb3MBM3/10._dp -        & 
&  9._dp*TrYx3adjYx3*Yb3adjYb3MBM3/5._dp - 9._dp*Yb3adjYb3Yb3adjYb3MBM3/50._dp -         & 
&  3._dp*Yb3adjYeYeadjYb3MBM3/5._dp - 9._dp*Yb3adjYw3Yw3adjYb3MBM3/10._dp

 
DMBM3 = oo16pi2*( betaMBM31 + oo16pi2 * betaMBM32 ) 

 
Else 
DMBM3 = oo16pi2* betaMBM31 
End If 
 
 
!-------------------- 
! TYu 
!-------------------- 
 
betaTYu1  = TYuadjYdYd + 5._dp*TYuadjYuYu + 26._dp*g1p2*MassB*Yu/15._dp +             & 
&  32._dp*g3p2*MassG*Yu/3._dp + 6._dp*g2p2*MassWB*Yu + 3._dp*TrTYb3adjYb3*Yu/5._dp +     & 
&  6._dp*TrTYuadjYu*Yu + 3._dp*TrTYw3adjYw3*Yu + 6._dp*TrTYx3adjYx3*Yu + 2._dp*YuadjYdTYd +& 
&  4._dp*YuadjYuTYu - 13._dp*g1p2*TYu/15._dp - 3._dp*g2p2*TYu - 16._dp*g3p2*TYu/3._dp +  & 
&  3._dp*TrYb3adjYb3*TYu/10._dp + 3._dp*TrYuadjYu*TYu + 3._dp*TrYw3adjYw3*TYu/2._dp +    & 
&  3._dp*TrYx3adjYx3*TYu

 
 
If (TwoLoopRGE) Then 
betaTYu2 = -2._dp*TYuadjYdTpYx3CYx3Yd + 2._dp*g1p2*TYuadjYdYd/5._dp - 3._dp*TrYdadjYd*TYuadjYdYd -& 
&  TrYeadjYe*TYuadjYdYd - 2._dp*TYuadjYdYdadjYdYd - 4._dp*TYuadjYdYdadjYuYu +            & 
&  12._dp*g2p2*TYuadjYuYu - 3._dp*TrYb3adjYb3*TYuadjYuYu/2._dp - 15._dp*TrYuadjYu*TYuadjYuYu -& 
&  15._dp*TrYw3adjYw3*TYuadjYuYu/2._dp - 15._dp*TrYx3adjYx3*TYuadjYuYu - 6._dp*TYuadjYuYuadjYuYu -& 
&  5486._dp*g1p4*MassB*Yu/225._dp - 2._dp*g1p2*g2p2*MassB*Yu - 272._dp*g1p2*g3p2*MassB*Yu/45._dp -& 
&  272._dp*g1p2*g3p2*MassG*Yu/45._dp - 16._dp*g2p2*g3p2*MassG*Yu + 64._dp*g3p4*MassG*Yu/9._dp -& 
&  2._dp*g1p2*g2p2*MassWB*Yu - 30._dp*g2p4*MassWB*Yu - 16._dp*g2p2*g3p2*MassWB*Yu -      & 
&  12._dp*TrCYdTpYdadjYx3TYx3*Yu - 12._dp*TrTpYx3CYx3TYdadjYd*Yu + 8._dp*g1p2*TrTYuadjYu*Yu/5._dp +& 
&  32._dp*g3p2*TrTYuadjYu*Yu + 12._dp*g2p2*TrTYw3adjYw3*Yu + 4._dp*g1p2*TrTYx3adjYx3*Yu +& 
&  32._dp*g3p2*TrTYx3adjYx3*Yu - 27._dp*TrYb3adjYb3TYb3adjYb3*Yu/25._dp - 3._dp*TrYb3adjYeTYeadjYb3*Yu/5._dp -& 
&  27._dp*TrYb3adjYw3TYw3adjYb3*Yu/10._dp - 6._dp*TrYdadjYuTYuadjYd*Yu - 3._dp*TrYeadjYb3TYb3adjYe*Yu/5._dp  
betaTYu2 =  betaTYu2- 3._dp*TrYeadjYw3TYw3adjYe*Yu - 6._dp*TrYuadjYdTYdadjYu*Yu - 8._dp*g1p2*MassB*TrYuadjYu*Yu/5._dp -& 
&  32._dp*g3p2*MassG*TrYuadjYu*Yu - 36._dp*TrYuadjYuTYuadjYu*Yu - 27._dp*TrYw3adjYb3TYb3adjYw3*Yu/10._dp -& 
&  3._dp*TrYw3adjYeTYeadjYw3*Yu - 12._dp*g2p2*MassWB*TrYw3adjYw3*Yu - 27._dp*TrYw3adjYw3TYw3adjYw3*Yu/2._dp -& 
&  4._dp*g1p2*MassB*TrYx3adjYx3*Yu - 32._dp*g3p2*MassG*TrYx3adjYx3*Yu - 36._dp*TrYx3adjYx3TYx3adjYx3*Yu -& 
&  4._dp*YuadjYdTpTYx3CYx3Yd - 4._dp*YuadjYdTpYx3CYx3TYd + 4._dp*g1p2*YuadjYdTYd/5._dp - & 
&  6._dp*TrYdadjYd*YuadjYdTYd - 2._dp*TrYeadjYe*YuadjYdTYd - 4._dp*YuadjYdTYdadjYdYd -   & 
&  4._dp*YuadjYdTYdadjYuYu - 4._dp*g1p2*MassB*YuadjYdYd/5._dp - 6._dp*TrTYdadjYd*YuadjYdYd -& 
&  2._dp*TrTYeadjYe*YuadjYdYd - 4._dp*YuadjYdYdadjYdTYd - 2._dp*YuadjYdYdadjYuTYu +      & 
&  6._dp*g1p2*YuadjYuTYu/5._dp + 6._dp*g2p2*YuadjYuTYu - 6._dp*TrYb3adjYb3*YuadjYuTYu/5._dp -& 
&  12._dp*TrYuadjYu*YuadjYuTYu - 6._dp*TrYw3adjYw3*YuadjYuTYu - 12._dp*TrYx3adjYx3*YuadjYuTYu -& 
&  8._dp*YuadjYuTYuadjYuYu - 4._dp*g1p2*MassB*YuadjYuYu/5._dp - 12._dp*g2p2*MassWB*YuadjYuYu  
betaTYu2 =  betaTYu2- 9._dp*TrTYb3adjYb3*YuadjYuYu/5._dp - 18._dp*TrTYuadjYu*YuadjYuYu - 9._dp*TrTYw3adjYw3*YuadjYuYu -& 
&  18._dp*TrTYx3adjYx3*YuadjYuYu - 6._dp*YuadjYuYuadjYuTYu - 64._dp*g3p4*MassG*Yu*NGHg3 -& 
&  24._dp*g2p4*MassWB*Yu*NGHw3 - 26._dp*g1p4*MassB*Yu*NGHx3/3._dp -& 
&  64._dp*g3p4*MassG*Yu*NGHx3/3._dp - 18._dp*g2p4*MassWB*Yu*NGHx3 -& 
&  26._dp*g1p4*MassB*Yu*NGHxb3/3._dp - 64._dp*g3p4*MassG*Yu*NGHxb3/3._dp -& 
&  18._dp*g2p4*MassWB*Yu*NGHxb3 + 2743._dp*g1p4*TYu/450._dp +             & 
&  g1p2*g2p2*TYu + 15._dp*g2p4*TYu/2._dp + 136._dp*g1p2*g3p2*TYu/45._dp + 8._dp*g2p2*g3p2*TYu -& 
&  16._dp*g3p4*TYu/9._dp - 6._dp*TrTpYx3CYx3YdadjYd*TYu - 27._dp*TrYb3adjYb3Yb3adjYb3*TYu/100._dp -& 
&  27._dp*TrYb3adjYw3Yw3adjYb3*TYu/40._dp - 3._dp*TrYdadjYuYuadjYd*TYu - 3._dp*TrYeadjYb3Yb3adjYe*TYu/10._dp -& 
&  3._dp*TrYeadjYw3Yw3adjYe*TYu/2._dp + 4._dp*g1p2*TrYuadjYu*TYu/5._dp + 16._dp*g3p2*TrYuadjYu*TYu -& 
&  9._dp*TrYuadjYuYuadjYu*TYu - 27._dp*TrYw3adjYb3Yb3adjYw3*TYu/40._dp + 6._dp*g2p2*TrYw3adjYw3*TYu  
betaTYu2 =  betaTYu2- 27._dp*TrYw3adjYw3Yw3adjYw3*TYu/8._dp + 2._dp*g1p2*TrYx3adjYx3*TYu + 16._dp*g3p2*TrYx3adjYx3*TYu -& 
&  9._dp*TrYx3adjYx3Yx3adjYx3*TYu + 16._dp*g3p4*NGHg3*TYu +               & 
&  6._dp*g2p4*NGHw3*TYu + 13._dp*g1p4*NGHx3*TYu/6._dp +    & 
&  9._dp*g2p4*NGHx3*TYu/2._dp + 16._dp*g3p4*NGHx3*TYu/3._dp +& 
&  13._dp*g1p4*NGHxb3*TYu/6._dp + 9._dp*g2p4*NGHxb3*TYu/2._dp +& 
&  16._dp*g3p4*NGHxb3*TYu/3._dp

 
DTYu = oo16pi2*( betaTYu1 + oo16pi2 * betaTYu2 ) 

 
Else 
DTYu = oo16pi2* betaTYu1 
End If 
 
 
!-------------------- 
! TYd 
!-------------------- 
 
betaTYd1  = 4._dp*TpTYx3CYx3Yd + 2._dp*TpYx3CYx3TYd + 5._dp*TYdadjYdYd +              & 
&  TYdadjYuYu + 14._dp*g1p2*MassB*Yd/15._dp + 32._dp*g3p2*MassG*Yd/3._dp +               & 
&  6._dp*g2p2*MassWB*Yd + 6._dp*TrTYdadjYd*Yd + 2._dp*TrTYeadjYe*Yd + 4._dp*YdadjYdTYd + & 
&  2._dp*YdadjYuTYu - 7._dp*g1p2*TYd/15._dp - 3._dp*g2p2*TYd - 16._dp*g3p2*TYd/3._dp +   & 
&  3._dp*TrYdadjYd*TYd + TrYeadjYe*TYd

 
 
If (TwoLoopRGE) Then 
betaTYd2 = -4._dp*TpTYx3CYx3TpYx3CYx3Yd + 4._dp*g1p2*TpTYx3CYx3Yd + 12._dp*g2p2*TpTYx3CYx3Yd -   & 
&  4._dp*TpYx3CYx3TpTYx3CYx3Yd - 2._dp*TpYx3CYx3TpYx3CYx3TYd + 2._dp*g1p2*TpYx3CYx3TYd + & 
&  6._dp*g2p2*TpYx3CYx3TYd - 4._dp*g1p2*MassB*TpYx3CYx3Yd - 12._dp*g2p2*MassWB*TpYx3CYx3Yd -& 
&  6._dp*TpYx3CYx3Yd*TrTYb3adjYb3/5._dp - 12._dp*TpYx3CYx3Yd*TrTYuadjYu - 6._dp*TpYx3CYx3Yd*TrTYw3adjYw3 -& 
&  12._dp*TpYx3CYx3Yd*TrTYx3adjYx3 - 6._dp*TpTYx3CYx3Yd*TrYb3adjYb3/5._dp -              & 
&  3._dp*TpYx3CYx3TYd*TrYb3adjYb3/5._dp - 12._dp*TpTYx3CYx3Yd*TrYuadjYu - 6._dp*TpYx3CYx3TYd*TrYuadjYu -& 
&  6._dp*TpTYx3CYx3Yd*TrYw3adjYw3 - 3._dp*TpYx3CYx3TYd*TrYw3adjYw3 - 12._dp*TpTYx3CYx3Yd*TrYx3adjYx3 -& 
&  6._dp*TpYx3CYx3TYd*TrYx3adjYx3 - 2._dp*TYdadjYdTpYx3CYx3Yd + 6._dp*g1p2*TYdadjYdYd/5._dp +& 
&  12._dp*g2p2*TYdadjYdYd - 15._dp*TrYdadjYd*TYdadjYdYd - 5._dp*TrYeadjYe*TYdadjYdYd -   & 
&  6._dp*TYdadjYdYdadjYdYd + 4._dp*g1p2*TYdadjYuYu/5._dp - 3._dp*TrYb3adjYb3*TYdadjYuYu/10._dp -& 
&  3._dp*TrYuadjYu*TYdadjYuYu - 3._dp*TrYw3adjYw3*TYdadjYuYu/2._dp - 3._dp*TrYx3adjYx3*TYdadjYuYu  
betaTYd2 =  betaTYd2- 4._dp*TYdadjYuYuadjYdYd - 2._dp*TYdadjYuYuadjYuYu - 574._dp*g1p4*MassB*Yd/45._dp -    & 
&  2._dp*g1p2*g2p2*MassB*Yd - 16._dp*g1p2*g3p2*MassB*Yd/9._dp - 16._dp*g1p2*g3p2*MassG*Yd/9._dp -& 
&  16._dp*g2p2*g3p2*MassG*Yd + 64._dp*g3p4*MassG*Yd/9._dp - 2._dp*g1p2*g2p2*MassWB*Yd -  & 
&  30._dp*g2p4*MassWB*Yd - 16._dp*g2p2*g3p2*MassWB*Yd - 12._dp*TrCYdTpYdadjYx3TYx3*Yd -  & 
&  4._dp*g1p2*TrTYdadjYd*Yd/5._dp + 32._dp*g3p2*TrTYdadjYd*Yd - 12._dp*TrTYdadjYdTpYx3CYx3*Yd +& 
&  12._dp*g1p2*TrTYeadjYe*Yd/5._dp - 3._dp*TrYb3adjYeTYeadjYb3*Yd/5._dp + 4._dp*g1p2*MassB*TrYdadjYd*Yd/5._dp -& 
&  32._dp*g3p2*MassG*TrYdadjYd*Yd - 36._dp*TrYdadjYdTYdadjYd*Yd - 6._dp*TrYdadjYuTYuadjYd*Yd -& 
&  3._dp*TrYeadjYb3TYb3adjYe*Yd/5._dp - 12._dp*g1p2*MassB*TrYeadjYe*Yd/5._dp -           & 
&  12._dp*TrYeadjYeTYeadjYe*Yd - 3._dp*TrYeadjYw3TYw3adjYe*Yd - 6._dp*TrYuadjYdTYdadjYu*Yd -& 
&  3._dp*TrYw3adjYeTYeadjYw3*Yd - 4._dp*YdadjYdTpTYx3CYx3Yd - 4._dp*YdadjYdTpYx3CYx3TYd +& 
&  6._dp*g1p2*YdadjYdTYd/5._dp + 6._dp*g2p2*YdadjYdTYd - 12._dp*TrYdadjYd*YdadjYdTYd  
betaTYd2 =  betaTYd2- 4._dp*TrYeadjYe*YdadjYdTYd - 8._dp*YdadjYdTYdadjYdYd - 8._dp*g1p2*MassB*YdadjYdYd/5._dp -& 
&  12._dp*g2p2*MassWB*YdadjYdYd - 18._dp*TrTYdadjYd*YdadjYdYd - 6._dp*TrTYeadjYe*YdadjYdYd -& 
&  6._dp*YdadjYdYdadjYdTYd + 8._dp*g1p2*YdadjYuTYu/5._dp - 3._dp*TrYb3adjYb3*YdadjYuTYu/5._dp -& 
&  6._dp*TrYuadjYu*YdadjYuTYu - 3._dp*TrYw3adjYw3*YdadjYuTYu - 6._dp*TrYx3adjYx3*YdadjYuTYu -& 
&  4._dp*YdadjYuTYuadjYdYd - 4._dp*YdadjYuTYuadjYuYu - 8._dp*g1p2*MassB*YdadjYuYu/5._dp -& 
&  3._dp*TrTYb3adjYb3*YdadjYuYu/5._dp - 6._dp*TrTYuadjYu*YdadjYuYu - 3._dp*TrTYw3adjYw3*YdadjYuYu -& 
&  6._dp*TrTYx3adjYx3*YdadjYuYu - 2._dp*YdadjYuYuadjYdTYd - 4._dp*YdadjYuYuadjYuTYu -    & 
&  64._dp*g3p4*MassG*Yd*NGHg3 - 24._dp*g2p4*MassWB*Yd*NGHw3 -& 
&  14._dp*g1p4*MassB*Yd*NGHx3/3._dp - 64._dp*g3p4*MassG*Yd*NGHx3/3._dp -& 
&  18._dp*g2p4*MassWB*Yd*NGHx3 - 14._dp*g1p4*MassB*Yd*NGHxb3/3._dp -& 
&  64._dp*g3p4*MassG*Yd*NGHxb3/3._dp - 18._dp*g2p4*MassWB*Yd*NGHxb3  
betaTYd2 =  betaTYd2+ 287._dp*g1p4*TYd/90._dp + g1p2*g2p2*TYd + 15._dp*g2p4*TYd/2._dp + 8._dp*g1p2*g3p2*TYd/9._dp +& 
&  8._dp*g2p2*g3p2*TYd - 16._dp*g3p4*TYd/9._dp - 3._dp*TrYb3adjYeYeadjYb3*TYd/10._dp -   & 
&  2._dp*g1p2*TrYdadjYd*TYd/5._dp + 16._dp*g3p2*TrYdadjYd*TYd - 6._dp*TrYdadjYdTpYx3CYx3*TYd -& 
&  9._dp*TrYdadjYdYdadjYd*TYd + 6._dp*g1p2*TrYeadjYe*TYd/5._dp - 3._dp*TrYeadjYeYeadjYe*TYd -& 
&  3._dp*TrYuadjYdYdadjYu*TYd - 3._dp*TrYw3adjYeYeadjYw3*TYd/2._dp + 16._dp*g3p4*NGHg3*TYd +& 
&  6._dp*g2p4*NGHw3*TYd + 7._dp*g1p4*NGHx3*TYd/6._dp +     & 
&  9._dp*g2p4*NGHx3*TYd/2._dp + 16._dp*g3p4*NGHx3*TYd/3._dp +& 
&  7._dp*g1p4*NGHxb3*TYd/6._dp + 9._dp*g2p4*NGHxb3*TYd/2._dp +& 
&  16._dp*g3p4*NGHxb3*TYd/3._dp

 
DTYd = oo16pi2*( betaTYd1 + oo16pi2 * betaTYd2 ) 

 
Else 
DTYd = oo16pi2* betaTYd1 
End If 
 
 
!-------------------- 
! TYe 
!-------------------- 
 
betaTYe1  = 3._dp*TYeadjYb3Yb3/10._dp + 5._dp*TYeadjYeYe + 3._dp*TYeadjYw3Yw3/2._dp + & 
&  18._dp*g1p2*MassB*Ye/5._dp + 6._dp*g2p2*MassWB*Ye + 6._dp*TrTYdadjYd*Ye +             & 
&  2._dp*TrTYeadjYe*Ye + 3._dp*YeadjYb3TYb3/5._dp + 4._dp*YeadjYeTYe + 3._dp*YeadjYw3TYw3 -& 
&  9._dp*g1p2*TYe/5._dp - 3._dp*g2p2*TYe + 3._dp*TrYdadjYd*TYe + TrYeadjYe*TYe

 
 
If (TwoLoopRGE) Then 
betaTYe2 = -9._dp*TrYb3adjYb3*TYeadjYb3Yb3/100._dp - 9._dp*TrYuadjYu*TYeadjYb3Yb3/10._dp -       & 
&  9._dp*TrYw3adjYw3*TYeadjYb3Yb3/20._dp - 9._dp*TrYx3adjYx3*TYeadjYb3Yb3/10._dp -       & 
&  9._dp*TYeadjYb3Yb3adjYb3Yb3/50._dp - 6._dp*TYeadjYb3Yb3adjYeYe/5._dp - 9._dp*TYeadjYb3Yb3adjYw3Yw3/40._dp -& 
&  6._dp*g1p2*TYeadjYeYe/5._dp + 12._dp*g2p2*TYeadjYeYe - 15._dp*TrYdadjYd*TYeadjYeYe -  & 
&  5._dp*TrYeadjYe*TYeadjYeYe - 6._dp*TYeadjYeYeadjYeYe + 6._dp*g2p2*TYeadjYw3Yw3 -      & 
&  9._dp*TrYb3adjYb3*TYeadjYw3Yw3/20._dp - 9._dp*TrYuadjYu*TYeadjYw3Yw3/2._dp -          & 
&  9._dp*TrYw3adjYw3*TYeadjYw3Yw3/4._dp - 9._dp*TrYx3adjYx3*TYeadjYw3Yw3/2._dp -         & 
&  9._dp*TYeadjYw3Yw3adjYb3Yb3/40._dp - 6._dp*TYeadjYw3Yw3adjYeYe - 9._dp*TYeadjYw3Yw3adjYw3Yw3/8._dp -& 
&  54._dp*g1p4*MassB*Ye - 18._dp*g1p2*g2p2*MassB*Ye/5._dp - 18._dp*g1p2*g2p2*MassWB*Ye/5._dp -& 
&  30._dp*g2p4*MassWB*Ye - 12._dp*TrCYdTpYdadjYx3TYx3*Ye - 4._dp*g1p2*TrTYdadjYd*Ye/5._dp +& 
&  32._dp*g3p2*TrTYdadjYd*Ye - 12._dp*TrTYdadjYdTpYx3CYx3*Ye + 12._dp*g1p2*TrTYeadjYe*Ye/5._dp  
betaTYe2 =  betaTYe2- 3._dp*TrYb3adjYeTYeadjYb3*Ye/5._dp + 4._dp*g1p2*MassB*TrYdadjYd*Ye/5._dp -            & 
&  32._dp*g3p2*MassG*TrYdadjYd*Ye - 36._dp*TrYdadjYdTYdadjYd*Ye - 6._dp*TrYdadjYuTYuadjYd*Ye -& 
&  3._dp*TrYeadjYb3TYb3adjYe*Ye/5._dp - 12._dp*g1p2*MassB*TrYeadjYe*Ye/5._dp -           & 
&  12._dp*TrYeadjYeTYeadjYe*Ye - 3._dp*TrYeadjYw3TYw3adjYe*Ye - 6._dp*TrYuadjYdTYdadjYu*Ye -& 
&  3._dp*TrYw3adjYeTYeadjYw3*Ye - 9._dp*TrYb3adjYb3*YeadjYb3TYb3/50._dp - 9._dp*TrYuadjYu*YeadjYb3TYb3/5._dp -& 
&  9._dp*TrYw3adjYw3*YeadjYb3TYb3/10._dp - 9._dp*TrYx3adjYx3*YeadjYb3TYb3/5._dp -        & 
&  9._dp*YeadjYb3TYb3adjYb3Yb3/25._dp - 6._dp*YeadjYb3TYb3adjYeYe/5._dp - 9._dp*YeadjYb3TYb3adjYw3Yw3/20._dp -& 
&  9._dp*TrTYb3adjYb3*YeadjYb3Yb3/50._dp - 9._dp*TrTYuadjYu*YeadjYb3Yb3/5._dp -          & 
&  9._dp*TrTYw3adjYw3*YeadjYb3Yb3/10._dp - 9._dp*TrTYx3adjYx3*YeadjYb3Yb3/5._dp -        & 
&  9._dp*YeadjYb3Yb3adjYb3TYb3/25._dp - 3._dp*YeadjYb3Yb3adjYeTYe/5._dp - 9._dp*YeadjYb3Yb3adjYw3TYw3/20._dp +& 
&  6._dp*g1p2*YeadjYeTYe/5._dp + 6._dp*g2p2*YeadjYeTYe - 12._dp*TrYdadjYd*YeadjYeTYe  
betaTYe2 =  betaTYe2- 4._dp*TrYeadjYe*YeadjYeTYe - 8._dp*YeadjYeTYeadjYeYe - 12._dp*g2p2*MassWB*YeadjYeYe - & 
&  18._dp*TrTYdadjYd*YeadjYeYe - 6._dp*TrTYeadjYe*YeadjYeYe - 6._dp*YeadjYeYeadjYeTYe +  & 
&  12._dp*g2p2*YeadjYw3TYw3 - 9._dp*TrYb3adjYb3*YeadjYw3TYw3/10._dp - 9._dp*TrYuadjYu*YeadjYw3TYw3 -& 
&  9._dp*TrYw3adjYw3*YeadjYw3TYw3/2._dp - 9._dp*TrYx3adjYx3*YeadjYw3TYw3 -               & 
&  9._dp*YeadjYw3TYw3adjYb3Yb3/20._dp - 6._dp*YeadjYw3TYw3adjYeYe - 9._dp*YeadjYw3TYw3adjYw3Yw3/4._dp -& 
&  12._dp*g2p2*MassWB*YeadjYw3Yw3 - 9._dp*TrTYb3adjYb3*YeadjYw3Yw3/10._dp -              & 
&  9._dp*TrTYuadjYu*YeadjYw3Yw3 - 9._dp*TrTYw3adjYw3*YeadjYw3Yw3/2._dp - 9._dp*TrTYx3adjYx3*YeadjYw3Yw3 -& 
&  9._dp*YeadjYw3Yw3adjYb3TYb3/20._dp - 3._dp*YeadjYw3Yw3adjYeTYe - 9._dp*YeadjYw3Yw3adjYw3TYw3/4._dp -& 
&  24._dp*g2p4*MassWB*Ye*NGHw3 - 18._dp*g1p4*MassB*Ye*NGHx3 -& 
&  18._dp*g2p4*MassWB*Ye*NGHx3 - 18._dp*g1p4*MassB*Ye*NGHxb3 -& 
&  18._dp*g2p4*MassWB*Ye*NGHxb3 + 27._dp*g1p4*TYe/2._dp + 9._dp*g1p2*g2p2*TYe/5._dp  
betaTYe2 =  betaTYe2+ 15._dp*g2p4*TYe/2._dp - 3._dp*TrYb3adjYeYeadjYb3*TYe/10._dp - 2._dp*g1p2*TrYdadjYd*TYe/5._dp +& 
&  16._dp*g3p2*TrYdadjYd*TYe - 6._dp*TrYdadjYdTpYx3CYx3*TYe - 9._dp*TrYdadjYdYdadjYd*TYe +& 
&  6._dp*g1p2*TrYeadjYe*TYe/5._dp - 3._dp*TrYeadjYeYeadjYe*TYe - 3._dp*TrYuadjYdYdadjYu*TYe -& 
&  3._dp*TrYw3adjYeYeadjYw3*TYe/2._dp + 6._dp*g2p4*NGHw3*TYe +            & 
&  9._dp*g1p4*NGHx3*TYe/2._dp + 9._dp*g2p4*NGHx3*TYe/2._dp +& 
&  9._dp*g1p4*NGHxb3*TYe/2._dp + 9._dp*g2p4*NGHxb3*TYe/2._dp

 
DTYe = oo16pi2*( betaTYe1 + oo16pi2 * betaTYe2 ) 

 
Else 
DTYe = oo16pi2* betaTYe1 
End If 
 
 
!-------------------- 
! TYb3 
!-------------------- 
 
betaTYb31  = 3._dp*TYb3adjYb3Yb3/2._dp + TYb3adjYeYe + 3._dp*TYb3adjYw3Yw3 +          & 
&  6._dp*g1p2*MassB*Yb3/5._dp + 6._dp*g2p2*MassWB*Yb3 + 3._dp*TrTYb3adjYb3*Yb3/5._dp +   & 
&  6._dp*TrTYuadjYu*Yb3 + 3._dp*TrTYw3adjYw3*Yb3 + 6._dp*TrTYx3adjYx3*Yb3 +              & 
&  6._dp*Yb3adjYb3TYb3/5._dp + 2._dp*Yb3adjYeTYe + 15._dp*Yb3adjYw3TYw3/4._dp -          & 
&  3._dp*g1p2*TYb3/5._dp - 3._dp*g2p2*TYb3 + 3._dp*TrYb3adjYb3*TYb3/10._dp +             & 
&  3._dp*TrYuadjYu*TYb3 + 3._dp*TrYw3adjYw3*TYb3/2._dp + 3._dp*TrYx3adjYx3*TYb3

 
 
If (TwoLoopRGE) Then 
betaTYb32 = 18._dp*g1p2*TYb3adjYb3Yb3/25._dp + 18._dp*g2p2*TYb3adjYb3Yb3/5._dp - 9._dp*TrYb3adjYb3*TYb3adjYb3Yb3/20._dp -& 
&  9._dp*TrYuadjYu*TYb3adjYb3Yb3/2._dp - 9._dp*TrYw3adjYw3*TYb3adjYb3Yb3/4._dp -         & 
&  9._dp*TrYx3adjYx3*TYb3adjYb3Yb3/2._dp - 27._dp*TYb3adjYb3Yb3adjYb3Yb3/50._dp -        & 
&  27._dp*TYb3adjYb3Yb3adjYw3Yw3/40._dp + 6._dp*g1p2*TYb3adjYeYe/5._dp - 3._dp*TrYdadjYd*TYb3adjYeYe -& 
&  TrYeadjYe*TYb3adjYeYe - 6._dp*TYb3adjYeYeadjYb3Yb3/5._dp - 2._dp*TYb3adjYeYeadjYeYe - & 
&  3._dp*TYb3adjYeYeadjYw3Yw3/2._dp + 9._dp*g1p2*TYb3adjYw3Yw3/10._dp + 21._dp*g2p2*TYb3adjYw3Yw3/2._dp -& 
&  9._dp*TrYb3adjYb3*TYb3adjYw3Yw3/10._dp - 9._dp*TrYuadjYu*TYb3adjYw3Yw3 -              & 
&  9._dp*TrYw3adjYw3*TYb3adjYw3Yw3/2._dp - 9._dp*TrYx3adjYx3*TYb3adjYw3Yw3 -             & 
&  81._dp*TYb3adjYw3Yw3adjYb3Yb3/40._dp - 27._dp*TYb3adjYw3Yw3adjYw3Yw3/8._dp -          & 
&  414._dp*g1p4*MassB*Yb3/25._dp - 18._dp*g1p2*g2p2*MassB*Yb3/5._dp - 18._dp*g1p2*g2p2*MassWB*Yb3/5._dp -& 
&  30._dp*g2p4*MassWB*Yb3 - 12._dp*TrCYdTpYdadjYx3TYx3*Yb3 - 12._dp*TrTpYx3CYx3TYdadjYd*Yb3  
betaTYb32 =  betaTYb32+ 8._dp*g1p2*TrTYuadjYu*Yb3/5._dp + 32._dp*g3p2*TrTYuadjYu*Yb3 + 12._dp*g2p2*TrTYw3adjYw3*Yb3 +& 
&  4._dp*g1p2*TrTYx3adjYx3*Yb3 + 32._dp*g3p2*TrTYx3adjYx3*Yb3 - 27._dp*TrYb3adjYb3TYb3adjYb3*Yb3/25._dp -& 
&  3._dp*TrYb3adjYeTYeadjYb3*Yb3/5._dp - 27._dp*TrYb3adjYw3TYw3adjYb3*Yb3/10._dp -       & 
&  6._dp*TrYdadjYuTYuadjYd*Yb3 - 3._dp*TrYeadjYb3TYb3adjYe*Yb3/5._dp - 3._dp*TrYeadjYw3TYw3adjYe*Yb3 -& 
&  6._dp*TrYuadjYdTYdadjYu*Yb3 - 8._dp*g1p2*MassB*TrYuadjYu*Yb3/5._dp - 32._dp*g3p2*MassG*TrYuadjYu*Yb3 -& 
&  36._dp*TrYuadjYuTYuadjYu*Yb3 - 27._dp*TrYw3adjYb3TYb3adjYw3*Yb3/10._dp -              & 
&  3._dp*TrYw3adjYeTYeadjYw3*Yb3 - 12._dp*g2p2*MassWB*TrYw3adjYw3*Yb3 - 27._dp*TrYw3adjYw3TYw3adjYw3*Yb3/2._dp -& 
&  4._dp*g1p2*MassB*TrYx3adjYx3*Yb3 - 32._dp*g3p2*MassG*TrYx3adjYx3*Yb3 - 36._dp*TrYx3adjYx3TYx3adjYx3*Yb3 +& 
&  9._dp*g1p2*Yb3adjYb3TYb3/25._dp + 9._dp*g2p2*Yb3adjYb3TYb3/5._dp - 9._dp*TrYb3adjYb3*Yb3adjYb3TYb3/25._dp -& 
&  18._dp*TrYuadjYu*Yb3adjYb3TYb3/5._dp - 9._dp*TrYw3adjYw3*Yb3adjYb3TYb3/5._dp -        & 
&  18._dp*TrYx3adjYx3*Yb3adjYb3TYb3/5._dp - 18._dp*Yb3adjYb3TYb3adjYb3Yb3/25._dp  
betaTYb32 =  betaTYb32- 9._dp*Yb3adjYb3TYb3adjYw3Yw3/10._dp - 18._dp*g1p2*MassB*Yb3adjYb3Yb3/25._dp -         & 
&  18._dp*g2p2*MassWB*Yb3adjYb3Yb3/5._dp - 27._dp*TrTYb3adjYb3*Yb3adjYb3Yb3/50._dp -     & 
&  27._dp*TrTYuadjYu*Yb3adjYb3Yb3/5._dp - 27._dp*TrTYw3adjYw3*Yb3adjYb3Yb3/10._dp -      & 
&  27._dp*TrTYx3adjYx3*Yb3adjYb3Yb3/5._dp - 27._dp*Yb3adjYb3Yb3adjYb3TYb3/50._dp -       & 
&  27._dp*Yb3adjYb3Yb3adjYw3TYw3/40._dp + 12._dp*g1p2*Yb3adjYeTYe/5._dp - 6._dp*TrYdadjYd*Yb3adjYeTYe -& 
&  2._dp*TrYeadjYe*Yb3adjYeTYe - 6._dp*Yb3adjYeTYeadjYb3Yb3/5._dp - 4._dp*Yb3adjYeTYeadjYeYe -& 
&  3._dp*Yb3adjYeTYeadjYw3Yw3/2._dp - 12._dp*g1p2*MassB*Yb3adjYeYe/5._dp -               & 
&  6._dp*TrTYdadjYd*Yb3adjYeYe - 2._dp*TrTYeadjYe*Yb3adjYeYe - 3._dp*Yb3adjYeYeadjYb3TYb3/5._dp -& 
&  4._dp*Yb3adjYeYeadjYeTYe - 3._dp*Yb3adjYeYeadjYw3TYw3/4._dp + 9._dp*g1p2*Yb3adjYw3TYw3/20._dp +& 
&  57._dp*g2p2*Yb3adjYw3TYw3/4._dp - 9._dp*TrYb3adjYb3*Yb3adjYw3TYw3/8._dp -             & 
&  45._dp*TrYuadjYu*Yb3adjYw3TYw3/4._dp - 45._dp*TrYw3adjYw3*Yb3adjYw3TYw3/8._dp  
betaTYb32 =  betaTYb32- 45._dp*TrYx3adjYx3*Yb3adjYw3TYw3/4._dp - 9._dp*Yb3adjYw3TYw3adjYb3Yb3/4._dp -         & 
&  9._dp*Yb3adjYw3TYw3adjYw3Yw3/2._dp - 9._dp*g1p2*MassB*Yb3adjYw3Yw3/10._dp -           & 
&  33._dp*g2p2*MassWB*Yb3adjYw3Yw3/2._dp - 27._dp*TrTYb3adjYb3*Yb3adjYw3Yw3/20._dp -     & 
&  27._dp*TrTYuadjYu*Yb3adjYw3Yw3/2._dp - 27._dp*TrTYw3adjYw3*Yb3adjYw3Yw3/4._dp -       & 
&  27._dp*TrTYx3adjYx3*Yb3adjYw3Yw3/2._dp - 27._dp*Yb3adjYw3Yw3adjYb3TYb3/20._dp -       & 
&  27._dp*Yb3adjYw3Yw3adjYw3TYw3/8._dp - 24._dp*g2p4*MassWB*Yb3*NGHw3 -   & 
&  6._dp*g1p4*MassB*Yb3*NGHx3 - 18._dp*g2p4*MassWB*Yb3*NGHx3 -& 
&  6._dp*g1p4*MassB*Yb3*NGHxb3 - 18._dp*g2p4*MassWB*Yb3*NGHxb3 +& 
&  207._dp*g1p4*TYb3/50._dp + 9._dp*g1p2*g2p2*TYb3/5._dp + 15._dp*g2p4*TYb3/2._dp -      & 
&  6._dp*TrTpYx3CYx3YdadjYd*TYb3 - 27._dp*TrYb3adjYb3Yb3adjYb3*TYb3/100._dp -            & 
&  27._dp*TrYb3adjYw3Yw3adjYb3*TYb3/40._dp - 3._dp*TrYdadjYuYuadjYd*TYb3  
betaTYb32 =  betaTYb32- 3._dp*TrYeadjYb3Yb3adjYe*TYb3/10._dp - 3._dp*TrYeadjYw3Yw3adjYe*TYb3/2._dp +          & 
&  4._dp*g1p2*TrYuadjYu*TYb3/5._dp + 16._dp*g3p2*TrYuadjYu*TYb3 - 9._dp*TrYuadjYuYuadjYu*TYb3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*TYb3/40._dp + 6._dp*g2p2*TrYw3adjYw3*TYb3 -               & 
&  27._dp*TrYw3adjYw3Yw3adjYw3*TYb3/8._dp + 2._dp*g1p2*TrYx3adjYx3*TYb3 + 16._dp*g3p2*TrYx3adjYx3*TYb3 -& 
&  9._dp*TrYx3adjYx3Yx3adjYx3*TYb3 + 6._dp*g2p4*NGHw3*TYb3 +              & 
&  3._dp*g1p4*NGHx3*TYb3/2._dp + 9._dp*g2p4*NGHx3*TYb3/2._dp +& 
&  3._dp*g1p4*NGHxb3*TYb3/2._dp + 9._dp*g2p4*NGHxb3*TYb3/2._dp

 
DTYb3 = oo16pi2*( betaTYb31 + oo16pi2 * betaTYb32 ) 

 
Else 
DTYb3 = oo16pi2* betaTYb31 
End If 
 
 
!-------------------- 
! TYw3 
!-------------------- 
 
betaTYw31  = 9._dp*TYw3adjYb3Yb3/10._dp + TYw3adjYeYe + 3._dp*TYw3adjYw3Yw3 +         & 
&  6._dp*g1p2*MassB*Yw3/5._dp + 14._dp*g2p2*MassWB*Yw3 + 3._dp*TrTYb3adjYb3*Yw3/5._dp +  & 
&  6._dp*TrTYuadjYu*Yw3 + 3._dp*TrTYw3adjYw3*Yw3 + 6._dp*TrTYx3adjYx3*Yw3 +              & 
&  9._dp*Yw3adjYb3TYb3/10._dp + 2._dp*Yw3adjYeTYe + 15._dp*Yw3adjYw3TYw3/4._dp -         & 
&  3._dp*g1p2*TYw3/5._dp - 7._dp*g2p2*TYw3 + 3._dp*TrYb3adjYb3*TYw3/10._dp +             & 
&  3._dp*TrYuadjYu*TYw3 + 3._dp*TrYw3adjYw3*TYw3/2._dp + 3._dp*TrYx3adjYx3*TYw3

 
 
If (TwoLoopRGE) Then 
betaTYw32 = 9._dp*g1p2*TYw3adjYb3Yb3/25._dp - 3._dp*g2p2*TYw3adjYb3Yb3/5._dp - 27._dp*TrYb3adjYb3*TYw3adjYb3Yb3/100._dp -& 
&  27._dp*TrYuadjYu*TYw3adjYb3Yb3/10._dp - 27._dp*TrYw3adjYw3*TYw3adjYb3Yb3/20._dp -     & 
&  27._dp*TrYx3adjYx3*TYw3adjYb3Yb3/10._dp - 9._dp*TYw3adjYb3Yb3adjYb3Yb3/25._dp -       & 
&  27._dp*TYw3adjYb3Yb3adjYw3Yw3/40._dp + 6._dp*g1p2*TYw3adjYeYe/5._dp - 3._dp*TrYdadjYd*TYw3adjYeYe -& 
&  TrYeadjYe*TYw3adjYeYe - 3._dp*TYw3adjYeYeadjYb3Yb3/5._dp - 2._dp*TYw3adjYeYeadjYeYe - & 
&  3._dp*TYw3adjYeYeadjYw3Yw3/2._dp + 9._dp*g1p2*TYw3adjYw3Yw3/10._dp + 9._dp*g2p2*TYw3adjYw3Yw3/2._dp -& 
&  9._dp*TrYb3adjYb3*TYw3adjYw3Yw3/10._dp - 9._dp*TrYuadjYu*TYw3adjYw3Yw3 -              & 
&  9._dp*TrYw3adjYw3*TYw3adjYw3Yw3/2._dp - 9._dp*TrYx3adjYx3*TYw3adjYw3Yw3 -             & 
&  9._dp*TYw3adjYw3Yw3adjYb3Yb3/8._dp - 27._dp*TYw3adjYw3Yw3adjYw3Yw3/8._dp -            & 
&  414._dp*g1p4*MassB*Yw3/25._dp - 18._dp*g1p2*g2p2*MassB*Yw3/5._dp - 18._dp*g1p2*g2p2*MassWB*Yw3/5._dp -& 
&  110._dp*g2p4*MassWB*Yw3 - 12._dp*TrCYdTpYdadjYx3TYx3*Yw3 - 12._dp*TrTpYx3CYx3TYdadjYd*Yw3  
betaTYw32 =  betaTYw32+ 8._dp*g1p2*TrTYuadjYu*Yw3/5._dp + 32._dp*g3p2*TrTYuadjYu*Yw3 + 12._dp*g2p2*TrTYw3adjYw3*Yw3 +& 
&  4._dp*g1p2*TrTYx3adjYx3*Yw3 + 32._dp*g3p2*TrTYx3adjYx3*Yw3 - 27._dp*TrYb3adjYb3TYb3adjYb3*Yw3/25._dp -& 
&  3._dp*TrYb3adjYeTYeadjYb3*Yw3/5._dp - 27._dp*TrYb3adjYw3TYw3adjYb3*Yw3/10._dp -       & 
&  6._dp*TrYdadjYuTYuadjYd*Yw3 - 3._dp*TrYeadjYb3TYb3adjYe*Yw3/5._dp - 3._dp*TrYeadjYw3TYw3adjYe*Yw3 -& 
&  6._dp*TrYuadjYdTYdadjYu*Yw3 - 8._dp*g1p2*MassB*TrYuadjYu*Yw3/5._dp - 32._dp*g3p2*MassG*TrYuadjYu*Yw3 -& 
&  36._dp*TrYuadjYuTYuadjYu*Yw3 - 27._dp*TrYw3adjYb3TYb3adjYw3*Yw3/10._dp -              & 
&  3._dp*TrYw3adjYeTYeadjYw3*Yw3 - 12._dp*g2p2*MassWB*TrYw3adjYw3*Yw3 - 27._dp*TrYw3adjYw3TYw3adjYw3*Yw3/2._dp -& 
&  4._dp*g1p2*MassB*TrYx3adjYx3*Yw3 - 32._dp*g3p2*MassG*TrYx3adjYx3*Yw3 - 36._dp*TrYx3adjYx3TYx3adjYx3*Yw3 +& 
&  9._dp*g1p2*Yw3adjYb3TYb3/50._dp - 3._dp*g2p2*Yw3adjYb3TYb3/10._dp - 27._dp*TrYb3adjYb3*Yw3adjYb3TYb3/100._dp -& 
&  27._dp*TrYuadjYu*Yw3adjYb3TYb3/10._dp - 27._dp*TrYw3adjYw3*Yw3adjYb3TYb3/20._dp -     & 
&  27._dp*TrYx3adjYx3*Yw3adjYb3TYb3/10._dp - 27._dp*Yw3adjYb3TYb3adjYb3Yb3/50._dp  
betaTYw32 =  betaTYw32- 9._dp*Yw3adjYb3TYb3adjYw3Yw3/10._dp - 9._dp*g1p2*MassB*Yw3adjYb3Yb3/25._dp +          & 
&  3._dp*g2p2*MassWB*Yw3adjYb3Yb3/5._dp - 9._dp*TrTYb3adjYb3*Yw3adjYb3Yb3/25._dp -       & 
&  18._dp*TrTYuadjYu*Yw3adjYb3Yb3/5._dp - 9._dp*TrTYw3adjYw3*Yw3adjYb3Yb3/5._dp -        & 
&  18._dp*TrTYx3adjYx3*Yw3adjYb3Yb3/5._dp - 9._dp*Yw3adjYb3Yb3adjYb3TYb3/20._dp -        & 
&  27._dp*Yw3adjYb3Yb3adjYw3TYw3/40._dp + 12._dp*g1p2*Yw3adjYeTYe/5._dp - 6._dp*TrYdadjYd*Yw3adjYeTYe -& 
&  2._dp*TrYeadjYe*Yw3adjYeTYe - 3._dp*Yw3adjYeTYeadjYb3Yb3/5._dp - 4._dp*Yw3adjYeTYeadjYeYe -& 
&  3._dp*Yw3adjYeTYeadjYw3Yw3/2._dp - 12._dp*g1p2*MassB*Yw3adjYeYe/5._dp -               & 
&  6._dp*TrTYdadjYd*Yw3adjYeYe - 2._dp*TrTYeadjYe*Yw3adjYeYe - 3._dp*Yw3adjYeYeadjYb3TYb3/10._dp -& 
&  4._dp*Yw3adjYeYeadjYeTYe - 3._dp*Yw3adjYeYeadjYw3TYw3/4._dp + 9._dp*g1p2*Yw3adjYw3TYw3/20._dp +& 
&  45._dp*g2p2*Yw3adjYw3TYw3/4._dp - 9._dp*TrYb3adjYb3*Yw3adjYw3TYw3/8._dp -             & 
&  45._dp*TrYuadjYu*Yw3adjYw3TYw3/4._dp - 45._dp*TrYw3adjYw3*Yw3adjYw3TYw3/8._dp  
betaTYw32 =  betaTYw32- 45._dp*TrYx3adjYx3*Yw3adjYw3TYw3/4._dp - 27._dp*Yw3adjYw3TYw3adjYb3Yb3/20._dp -       & 
&  9._dp*Yw3adjYw3TYw3adjYw3Yw3/2._dp - 9._dp*g1p2*MassB*Yw3adjYw3Yw3/10._dp -           & 
&  21._dp*g2p2*MassWB*Yw3adjYw3Yw3/2._dp - 27._dp*TrTYb3adjYb3*Yw3adjYw3Yw3/20._dp -     & 
&  27._dp*TrTYuadjYu*Yw3adjYw3Yw3/2._dp - 27._dp*TrTYw3adjYw3*Yw3adjYw3Yw3/4._dp -       & 
&  27._dp*TrTYx3adjYx3*Yw3adjYw3Yw3/2._dp - 9._dp*Yw3adjYw3Yw3adjYb3TYb3/10._dp -        & 
&  27._dp*Yw3adjYw3Yw3adjYw3TYw3/8._dp - 56._dp*g2p4*MassWB*Yw3*NGHw3 -   & 
&  6._dp*g1p4*MassB*Yw3*NGHx3 - 42._dp*g2p4*MassWB*Yw3*NGHx3 -& 
&  6._dp*g1p4*MassB*Yw3*NGHxb3 - 42._dp*g2p4*MassWB*Yw3*NGHxb3 +& 
&  207._dp*g1p4*TYw3/50._dp + 9._dp*g1p2*g2p2*TYw3/5._dp + 55._dp*g2p4*TYw3/2._dp -      & 
&  6._dp*TrTpYx3CYx3YdadjYd*TYw3 - 27._dp*TrYb3adjYb3Yb3adjYb3*TYw3/100._dp -            & 
&  27._dp*TrYb3adjYw3Yw3adjYb3*TYw3/40._dp - 3._dp*TrYdadjYuYuadjYd*TYw3  
betaTYw32 =  betaTYw32- 3._dp*TrYeadjYb3Yb3adjYe*TYw3/10._dp - 3._dp*TrYeadjYw3Yw3adjYe*TYw3/2._dp +          & 
&  4._dp*g1p2*TrYuadjYu*TYw3/5._dp + 16._dp*g3p2*TrYuadjYu*TYw3 - 9._dp*TrYuadjYuYuadjYu*TYw3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*TYw3/40._dp + 6._dp*g2p2*TrYw3adjYw3*TYw3 -               & 
&  27._dp*TrYw3adjYw3Yw3adjYw3*TYw3/8._dp + 2._dp*g1p2*TrYx3adjYx3*TYw3 + 16._dp*g3p2*TrYx3adjYx3*TYw3 -& 
&  9._dp*TrYx3adjYx3Yx3adjYx3*TYw3 + 14._dp*g2p4*NGHw3*TYw3 +             & 
&  3._dp*g1p4*NGHx3*TYw3/2._dp + 21._dp*g2p4*NGHx3*TYw3/2._dp +& 
&  3._dp*g1p4*NGHxb3*TYw3/2._dp + 21._dp*g2p4*NGHxb3*TYw3/2._dp

 
DTYw3 = oo16pi2*( betaTYw31 + oo16pi2 * betaTYw32 ) 

 
Else 
DTYw3 = oo16pi2* betaTYw31 
End If 
 
 
!-------------------- 
! TYx3 
!-------------------- 
 
betaTYx31  = 4._dp*TYx3adjYx3Yx3 + 2._dp*TYx3CYdTpYd + 38._dp*g1p2*MassB*Yx3/15._dp + & 
&  32._dp*g3p2*MassG*Yx3/3._dp + 6._dp*g2p2*MassWB*Yx3 + 3._dp*TrTYb3adjYb3*Yx3/5._dp +  & 
&  6._dp*TrTYuadjYu*Yx3 + 3._dp*TrTYw3adjYw3*Yx3 + 6._dp*TrTYx3adjYx3*Yx3 +              & 
&  5._dp*Yx3adjYx3TYx3 + 4._dp*Yx3CYdTpTYd - 19._dp*g1p2*TYx3/15._dp - 3._dp*g2p2*TYx3 - & 
&  16._dp*g3p2*TYx3/3._dp + 3._dp*TrYb3adjYb3*TYx3/10._dp + 3._dp*TrYuadjYu*TYx3 +       & 
&  3._dp*TrYw3adjYw3*TYx3/2._dp + 3._dp*TrYx3adjYx3*TYx3

 
 
If (TwoLoopRGE) Then 
betaTYx32 = 6._dp*g1p2*TYx3adjYx3Yx3/5._dp + 6._dp*g2p2*TYx3adjYx3Yx3 - 6._dp*TrYb3adjYb3*TYx3adjYx3Yx3/5._dp -& 
&  12._dp*TrYuadjYu*TYx3adjYx3Yx3 - 6._dp*TrYw3adjYw3*TYx3adjYx3Yx3 - 12._dp*TrYx3adjYx3*TYx3adjYx3Yx3 -& 
&  6._dp*TYx3adjYx3Yx3adjYx3Yx3 + 2._dp*g1p2*TYx3CYdTpYd/5._dp + 6._dp*g2p2*TYx3CYdTpYd -& 
&  6._dp*TrYdadjYd*TYx3CYdTpYd - 2._dp*TrYeadjYe*TYx3CYdTpYd - 4._dp*TYx3CYdTpYdadjYx3Yx3 -& 
&  2._dp*TYx3CYdTpYdCYdTpYd - 2._dp*TYx3CYdTpYuCYuTpYd - 8246._dp*g1p4*MassB*Yx3/225._dp -& 
&  34._dp*g1p2*g2p2*MassB*Yx3/5._dp - 464._dp*g1p2*g3p2*MassB*Yx3/45._dp -               & 
&  464._dp*g1p2*g3p2*MassG*Yx3/45._dp - 16._dp*g2p2*g3p2*MassG*Yx3 + 64._dp*g3p4*MassG*Yx3/9._dp -& 
&  34._dp*g1p2*g2p2*MassWB*Yx3/5._dp - 30._dp*g2p4*MassWB*Yx3 - 16._dp*g2p2*g3p2*MassWB*Yx3 -& 
&  12._dp*TrCYdTpYdadjYx3TYx3*Yx3 - 12._dp*TrTpYx3CYx3TYdadjYd*Yx3 + 8._dp*g1p2*TrTYuadjYu*Yx3/5._dp +& 
&  32._dp*g3p2*TrTYuadjYu*Yx3 + 12._dp*g2p2*TrTYw3adjYw3*Yx3 + 4._dp*g1p2*TrTYx3adjYx3*Yx3 +& 
&  32._dp*g3p2*TrTYx3adjYx3*Yx3 - 27._dp*TrYb3adjYb3TYb3adjYb3*Yx3/25._dp  
betaTYx32 =  betaTYx32- 3._dp*TrYb3adjYeTYeadjYb3*Yx3/5._dp - 27._dp*TrYb3adjYw3TYw3adjYb3*Yx3/10._dp -       & 
&  6._dp*TrYdadjYuTYuadjYd*Yx3 - 3._dp*TrYeadjYb3TYb3adjYe*Yx3/5._dp - 3._dp*TrYeadjYw3TYw3adjYe*Yx3 -& 
&  6._dp*TrYuadjYdTYdadjYu*Yx3 - 8._dp*g1p2*MassB*TrYuadjYu*Yx3/5._dp - 32._dp*g3p2*MassG*TrYuadjYu*Yx3 -& 
&  36._dp*TrYuadjYuTYuadjYu*Yx3 - 27._dp*TrYw3adjYb3TYb3adjYw3*Yx3/10._dp -              & 
&  3._dp*TrYw3adjYeTYeadjYw3*Yx3 - 12._dp*g2p2*MassWB*TrYw3adjYw3*Yx3 - 27._dp*TrYw3adjYw3TYw3adjYw3*Yx3/2._dp -& 
&  4._dp*g1p2*MassB*TrYx3adjYx3*Yx3 - 32._dp*g3p2*MassG*TrYx3adjYx3*Yx3 - 36._dp*TrYx3adjYx3TYx3adjYx3*Yx3 +& 
&  18._dp*g1p2*Yx3adjYx3TYx3/5._dp + 12._dp*g2p2*Yx3adjYx3TYx3 - 3._dp*TrYb3adjYb3*Yx3adjYx3TYx3/2._dp -& 
&  15._dp*TrYuadjYu*Yx3adjYx3TYx3 - 15._dp*TrYw3adjYw3*Yx3adjYx3TYx3/2._dp -             & 
&  15._dp*TrYx3adjYx3*Yx3adjYx3TYx3 - 8._dp*Yx3adjYx3TYx3adjYx3Yx3 - 16._dp*g1p2*MassB*Yx3adjYx3Yx3/5._dp -& 
&  12._dp*g2p2*MassWB*Yx3adjYx3Yx3 - 9._dp*TrTYb3adjYb3*Yx3adjYx3Yx3/5._dp -             & 
&  18._dp*TrTYuadjYu*Yx3adjYx3Yx3 - 9._dp*TrTYw3adjYw3*Yx3adjYx3Yx3 - 18._dp*TrTYx3adjYx3*Yx3adjYx3Yx3  
betaTYx32 =  betaTYx32- 6._dp*Yx3adjYx3Yx3adjYx3TYx3 + 4._dp*g1p2*Yx3CYdTpTYd/5._dp + 12._dp*g2p2*Yx3CYdTpTYd -& 
&  12._dp*TrYdadjYd*Yx3CYdTpTYd - 4._dp*TrYeadjYe*Yx3CYdTpTYd - 4._dp*Yx3CYdTpTYdCYdTpYd -& 
&  4._dp*Yx3CYdTpTYuCYuTpYd - 4._dp*g1p2*MassB*Yx3CYdTpYd/5._dp - 12._dp*g2p2*MassWB*Yx3CYdTpYd -& 
&  12._dp*TrTYdadjYd*Yx3CYdTpYd - 4._dp*TrTYeadjYe*Yx3CYdTpYd - 2._dp*Yx3CYdTpYdadjYx3TYx3 -& 
&  4._dp*Yx3CYdTpYdCYdTpTYd - 4._dp*Yx3CYdTpYuCYuTpTYd - 4._dp*Yx3TYdadjYdadjYx3Yx3 -    & 
&  64._dp*g3p4*MassG*Yx3*NGHg3 - 24._dp*g2p4*MassWB*Yx3*NGHw3 -& 
&  38._dp*g1p4*MassB*Yx3*NGHx3/3._dp - 64._dp*g3p4*MassG*Yx3*NGHx3/3._dp -& 
&  18._dp*g2p4*MassWB*Yx3*NGHx3 - 38._dp*g1p4*MassB*Yx3*NGHxb3/3._dp -& 
&  64._dp*g3p4*MassG*Yx3*NGHxb3/3._dp - 18._dp*g2p4*MassWB*Yx3*NGHxb3 +& 
&  4123._dp*g1p4*TYx3/450._dp + 17._dp*g1p2*g2p2*TYx3/5._dp + 15._dp*g2p4*TYx3/2._dp +   & 
&  232._dp*g1p2*g3p2*TYx3/45._dp + 8._dp*g2p2*g3p2*TYx3 - 16._dp*g3p4*TYx3/9._dp  
betaTYx32 =  betaTYx32- 6._dp*TrTpYx3CYx3YdadjYd*TYx3 - 27._dp*TrYb3adjYb3Yb3adjYb3*TYx3/100._dp -            & 
&  27._dp*TrYb3adjYw3Yw3adjYb3*TYx3/40._dp - 3._dp*TrYdadjYuYuadjYd*TYx3 -               & 
&  3._dp*TrYeadjYb3Yb3adjYe*TYx3/10._dp - 3._dp*TrYeadjYw3Yw3adjYe*TYx3/2._dp +          & 
&  4._dp*g1p2*TrYuadjYu*TYx3/5._dp + 16._dp*g3p2*TrYuadjYu*TYx3 - 9._dp*TrYuadjYuYuadjYu*TYx3 -& 
&  27._dp*TrYw3adjYb3Yb3adjYw3*TYx3/40._dp + 6._dp*g2p2*TrYw3adjYw3*TYx3 -               & 
&  27._dp*TrYw3adjYw3Yw3adjYw3*TYx3/8._dp + 2._dp*g1p2*TrYx3adjYx3*TYx3 + 16._dp*g3p2*TrYx3adjYx3*TYx3 -& 
&  9._dp*TrYx3adjYx3Yx3adjYx3*TYx3 + 16._dp*g3p4*NGHg3*TYx3 +             & 
&  6._dp*g2p4*NGHw3*TYx3 + 19._dp*g1p4*NGHx3*TYx3/6._dp +  & 
&  9._dp*g2p4*NGHx3*TYx3/2._dp + 16._dp*g3p4*NGHx3*TYx3/3._dp +& 
&  19._dp*g1p4*NGHxb3*TYx3/6._dp + 9._dp*g2p4*NGHxb3*TYx3/2._dp +& 
&  16._dp*g3p4*NGHxb3*TYx3/3._dp

 
DTYx3 = oo16pi2*( betaTYx31 + oo16pi2 * betaTYx32 ) 

 
Else 
DTYx3 = oo16pi2* betaTYx31 
End If 
 
 
!-------------------- 
! Bmue 
!-------------------- 
 
betaBmue1  = 6._dp*g1p2*MassB*mue/5._dp + 6._dp*g2p2*MassWB*mue + 3._dp*TrTYb3adjYb3*mue/5._dp +& 
&  6._dp*TrTYdadjYd*mue + 2._dp*TrTYeadjYe*mue + 6._dp*TrTYuadjYu*mue + 3._dp*TrTYw3adjYw3*mue +& 
&  6._dp*TrTYx3adjYx3*mue - 3._dp*g1p2*Bmue/5._dp - 3._dp*g2p2*Bmue + 3._dp*TrYb3adjYb3*Bmue/10._dp +& 
&  3._dp*TrYdadjYd*Bmue + TrYeadjYe*Bmue + 3._dp*TrYuadjYu*Bmue + 3._dp*TrYw3adjYw3*Bmue/2._dp +& 
&  3._dp*TrYx3adjYx3*Bmue

 
 
If (TwoLoopRGE) Then 
betaBmue2 = -414._dp*g1p4*MassB*mue/25._dp - 18._dp*g1p2*g2p2*MassB*mue/5._dp - 18._dp*g1p2*g2p2*MassWB*mue/5._dp -& 
&  30._dp*g2p4*MassWB*mue - 24._dp*TrCYdTpYdadjYx3TYx3*mue - 12._dp*TrTpYx3CYx3TYdadjYd*mue -& 
&  4._dp*g1p2*TrTYdadjYd*mue/5._dp + 32._dp*g3p2*TrTYdadjYd*mue - 12._dp*TrTYdadjYdTpYx3CYx3*mue +& 
&  12._dp*g1p2*TrTYeadjYe*mue/5._dp + 8._dp*g1p2*TrTYuadjYu*mue/5._dp + 32._dp*g3p2*TrTYuadjYu*mue +& 
&  12._dp*g2p2*TrTYw3adjYw3*mue + 4._dp*g1p2*TrTYx3adjYx3*mue + 32._dp*g3p2*TrTYx3adjYx3*mue -& 
&  27._dp*TrYb3adjYb3TYb3adjYb3*mue/25._dp - 6._dp*TrYb3adjYeTYeadjYb3*mue/5._dp -       & 
&  27._dp*TrYb3adjYw3TYw3adjYb3*mue/10._dp + 4._dp*g1p2*MassB*TrYdadjYd*mue/5._dp -      & 
&  32._dp*g3p2*MassG*TrYdadjYd*mue - 36._dp*TrYdadjYdTYdadjYd*mue - 12._dp*TrYdadjYuTYuadjYd*mue -& 
&  6._dp*TrYeadjYb3TYb3adjYe*mue/5._dp - 12._dp*g1p2*MassB*TrYeadjYe*mue/5._dp -         & 
&  12._dp*TrYeadjYeTYeadjYe*mue - 6._dp*TrYeadjYw3TYw3adjYe*mue - 12._dp*TrYuadjYdTYdadjYu*mue -& 
&  8._dp*g1p2*MassB*TrYuadjYu*mue/5._dp - 32._dp*g3p2*MassG*TrYuadjYu*mue  
betaBmue2 =  betaBmue2- 36._dp*TrYuadjYuTYuadjYu*mue - 27._dp*TrYw3adjYb3TYb3adjYw3*mue/10._dp -              & 
&  6._dp*TrYw3adjYeTYeadjYw3*mue - 12._dp*g2p2*MassWB*TrYw3adjYw3*mue - 27._dp*TrYw3adjYw3TYw3adjYw3*mue/2._dp -& 
&  4._dp*g1p2*MassB*TrYx3adjYx3*mue - 32._dp*g3p2*MassG*TrYx3adjYx3*mue - 36._dp*TrYx3adjYx3TYx3adjYx3*mue +& 
&  207._dp*g1p4*Bmue/50._dp + 9._dp*g1p2*g2p2*Bmue/5._dp + 15._dp*g2p4*Bmue/2._dp -      & 
&  6._dp*TrTpYx3CYx3YdadjYd*Bmue - 27._dp*TrYb3adjYb3Yb3adjYb3*Bmue/100._dp -            & 
&  3._dp*TrYb3adjYeYeadjYb3*Bmue/10._dp - 27._dp*TrYb3adjYw3Yw3adjYb3*Bmue/40._dp -      & 
&  2._dp*g1p2*TrYdadjYd*Bmue/5._dp + 16._dp*g3p2*TrYdadjYd*Bmue - 6._dp*TrYdadjYdTpYx3CYx3*Bmue -& 
&  9._dp*TrYdadjYdYdadjYd*Bmue - 3._dp*TrYdadjYuYuadjYd*Bmue - 3._dp*TrYeadjYb3Yb3adjYe*Bmue/10._dp +& 
&  6._dp*g1p2*TrYeadjYe*Bmue/5._dp - 3._dp*TrYeadjYeYeadjYe*Bmue - 3._dp*TrYeadjYw3Yw3adjYe*Bmue/2._dp -& 
&  3._dp*TrYuadjYdYdadjYu*Bmue + 4._dp*g1p2*TrYuadjYu*Bmue/5._dp + 16._dp*g3p2*TrYuadjYu*Bmue -& 
&  9._dp*TrYuadjYuYuadjYu*Bmue - 27._dp*TrYw3adjYb3Yb3adjYw3*Bmue/40._dp  

betaBmue2 =  betaBmue2- 3._dp*TrYw3adjYeYeadjYw3*Bmue/2._dp + 6._dp*g2p2*TrYw3adjYw3*Bmue - 27._dp*TrYw3adjYw3Yw3adjYw3*Bmue/8._dp+&
&  2._dp*g1p2*TrYx3adjYx3*Bmue + 16._dp*g3p2*TrYx3adjYx3*Bmue - 9._dp*TrYx3adjYx3Yx3adjYx3*Bmue -& 
&  24._dp*g2p4*MassWB*mue*NGHw3 + 6._dp*g2p4*Bmue*NGHw3 -  & 
&  6._dp*g1p4*MassB*mue*NGHx3 - 18._dp*g2p4*MassWB*mue*NGHx3 +& 
&  3._dp*g1p4*Bmue*NGHx3/2._dp + 9._dp*g2p4*Bmue*NGHx3/2._dp -& 
&  6._dp*g1p4*MassB*mue*NGHxb3 - 18._dp*g2p4*MassWB*mue*NGHxb3 +& 
&  3._dp*g1p4*Bmue*NGHxb3/2._dp + 9._dp*g2p4*Bmue*NGHxb3/2._dp

 
DBmue = oo16pi2*( betaBmue1 + oo16pi2 * betaBmue2 ) 

 
Else 
DBmue = oo16pi2* betaBmue1 
End If 
 
 
!-------------------- 
! BMXM3 
!-------------------- 
 
betaBMXM31  = BMXM3CYx3TpYx3 + 10._dp*g1p2*MassB*MXM3/3._dp + 32._dp*g3p2*MassG*MXM3/3._dp +& 
&  6._dp*g2p2*MassWB*MXM3 + 2._dp*MXM3CYx3TpTYx3 - 5._dp*g1p2*BMXM3/3._dp -              & 
&  3._dp*g2p2*BMXM3 - 16._dp*g3p2*BMXM3/3._dp

 
 
If (TwoLoopRGE) Then 
betaBMXM32 = -2._dp*BMXM3CYx3TpYx3CYx3TpYx3 - 2._dp*BMXM3CYx3YdadjYdTpYx3 - 2._dp*BMXM3CYx3TpYx3*g1p2/5._dp -& 
&  446._dp*g1p4*MassB*MXM3/9._dp - 10._dp*g1p2*g2p2*MassB*MXM3 - 160._dp*g1p2*g3p2*MassB*MXM3/9._dp -& 
&  160._dp*g1p2*g3p2*MassG*MXM3/9._dp - 32._dp*g2p2*g3p2*MassG*MXM3 + 64._dp*g3p4*MassG*MXM3/9._dp -& 
&  10._dp*g1p2*g2p2*MassWB*MXM3 - 30._dp*g2p4*MassWB*MXM3 - 32._dp*g2p2*g3p2*MassWB*MXM3 -& 
&  4._dp*g1p2*MXM3CYx3TpTYx3/5._dp - 4._dp*MXM3CYx3TpTYx3CYx3TpYx3 + 4._dp*g1p2*MassB*MXM3CYx3TpYx3/5._dp -& 
&  4._dp*MXM3CYx3TpYx3CYx3TpTYx3 - 4._dp*MXM3CYx3TYdadjYdTpYx3 - 4._dp*MXM3CYx3YdadjYdTpTYx3 -& 
&  3._dp*MXM3CYx3TpYx3*TrTYb3adjYb3/5._dp - 6._dp*MXM3CYx3TpYx3*TrTYuadjYu -             & 
&  3._dp*MXM3CYx3TpYx3*TrTYw3adjYw3 - 6._dp*MXM3CYx3TpYx3*TrTYx3adjYx3 - 3._dp*BMXM3CYx3TpYx3*TrYb3adjYb3/10._dp -& 
&  3._dp*MXM3CYx3TpTYx3*TrYb3adjYb3/5._dp - 3._dp*BMXM3CYx3TpYx3*TrYuadjYu -             & 
&  6._dp*MXM3CYx3TpTYx3*TrYuadjYu - 3._dp*BMXM3CYx3TpYx3*TrYw3adjYw3/2._dp -             & 
&  3._dp*MXM3CYx3TpTYx3*TrYw3adjYw3 - 3._dp*BMXM3CYx3TpYx3*TrYx3adjYx3 - 6._dp*MXM3CYx3TpTYx3*TrYx3adjYx3  
betaBMXM32 =  betaBMXM32+ 223._dp*g1p4*BMXM3/18._dp + 5._dp*g1p2*g2p2*BMXM3 + 15._dp*g2p4*BMXM3/2._dp +         & 
&  80._dp*g1p2*g3p2*BMXM3/9._dp + 16._dp*g2p2*g3p2*BMXM3 - 16._dp*g3p4*BMXM3/9._dp -     & 
&  64._dp*g3p4*MassG*MXM3*NGHg3 + 16._dp*g3p4*BMXM3*NGHg3 -& 
&  24._dp*g2p4*MassWB*MXM3*NGHw3 + 6._dp*g2p4*BMXM3*NGHw3 -& 
&  50._dp*g1p4*MassB*MXM3*NGHx3/3._dp - 64._dp*g3p4*MassG*MXM3*NGHx3/3._dp -& 
&  18._dp*g2p4*MassWB*MXM3*NGHx3 + 25._dp*g1p4*BMXM3*NGHx3/6._dp +& 
&  9._dp*g2p4*BMXM3*NGHx3/2._dp + 16._dp*g3p4*BMXM3*NGHx3/3._dp -& 
&  50._dp*g1p4*MassB*MXM3*NGHxb3/3._dp - 64._dp*g3p4*MassG*MXM3*NGHxb3/3._dp -& 
&  18._dp*g2p4*MassWB*MXM3*NGHxb3 + 25._dp*g1p4*BMXM3*NGHxb3/6._dp +& 
&  9._dp*g2p4*BMXM3*NGHxb3/2._dp + 16._dp*g3p4*BMXM3*NGHxb3/3._dp

 
DBMXM3 = oo16pi2*( betaBMXM31 + oo16pi2 * betaBMXM32 ) 

 
Else 
DBMXM3 = oo16pi2* betaBMXM31 
End If 
 
 
!-------------------- 
! BMWM3 
!-------------------- 
 
betaBMWM31  = BMWM3CYw3TpYw3 + 16._dp*g2p2*MassWB*MWM3 + 2._dp*MWM3CYw3TpTYw3 +       & 
&  2._dp*TYw3adjYw3MWM3 + Yw3adjYw3BMWM3 - 8._dp*g2p2*BMWM3

 
 
If (TwoLoopRGE) Then 
betaBMWM32 = -3._dp*BMWM3CYw3TpYb3CYb3TpYw3/10._dp - BMWM3CYw3TpYeCYeTpYw3 - 3._dp*BMWM3CYw3TpYw3CYw3TpYw3/2._dp +& 
&  3._dp*BMWM3CYw3TpYw3*g1p2/5._dp - BMWM3CYw3TpYw3*g2p2 - 160._dp*g2p4*MassWB*MWM3 -    & 
&  3._dp*MWM3CYw3TpTYb3CYb3TpYw3/5._dp - 2._dp*MWM3CYw3TpTYeCYeTpYw3 + 6._dp*g1p2*MWM3CYw3TpTYw3/5._dp -& 
&  2._dp*g2p2*MWM3CYw3TpTYw3 - 3._dp*MWM3CYw3TpTYw3CYw3TpYw3 - 3._dp*MWM3CYw3TpYb3CYb3TpTYw3/5._dp -& 
&  2._dp*MWM3CYw3TpYeCYeTpTYw3 - 6._dp*g1p2*MassB*MWM3CYw3TpYw3/5._dp + 2._dp*g2p2*MassWB*MWM3CYw3TpYw3 -& 
&  3._dp*MWM3CYw3TpYw3CYw3TpTYw3 - 3._dp*MWM3CYw3TpYw3*TrTYb3adjYb3/5._dp -              & 
&  6._dp*MWM3CYw3TpYw3*TrTYuadjYu - 3._dp*MWM3CYw3TpYw3*TrTYw3adjYw3 - 6._dp*MWM3CYw3TpYw3*TrTYx3adjYx3 -& 
&  3._dp*BMWM3CYw3TpYw3*TrYb3adjYb3/10._dp - 3._dp*MWM3CYw3TpTYw3*TrYb3adjYb3/5._dp -    & 
&  3._dp*BMWM3CYw3TpYw3*TrYuadjYu - 6._dp*MWM3CYw3TpTYw3*TrYuadjYu - 3._dp*BMWM3CYw3TpYw3*TrYw3adjYw3/2._dp -& 
&  3._dp*MWM3CYw3TpTYw3*TrYw3adjYw3 - 3._dp*BMWM3CYw3TpYw3*TrYx3adjYx3 - 6._dp*MWM3CYw3TpTYw3*TrYx3adjYx3 -& 
&  3._dp*TYw3adjYb3Yb3adjYw3MWM3/5._dp - 2._dp*TYw3adjYeYeadjYw3MWM3 + 6._dp*g1p2*TYw3adjYw3MWM3/5._dp  
betaBMWM32 =  betaBMWM32- 2._dp*g2p2*TYw3adjYw3MWM3 - 3._dp*TrYb3adjYb3*TYw3adjYw3MWM3/5._dp - 6._dp*TrYuadjYu*TYw3adjYw3MWM3 -& 
&  3._dp*TrYw3adjYw3*TYw3adjYw3MWM3 - 6._dp*TrYx3adjYx3*TYw3adjYw3MWM3 - 3._dp*TYw3adjYw3Yw3adjYw3MWM3 -& 
&  3._dp*Yw3adjYb3TYb3adjYw3MWM3/5._dp - 3._dp*Yw3adjYb3Yb3adjYw3BMWM3/10._dp -          & 
&  2._dp*Yw3adjYeTYeadjYw3MWM3 - Yw3adjYeYeadjYw3BMWM3 + 3._dp*g1p2*Yw3adjYw3BMWM3/5._dp -& 
&  g2p2*Yw3adjYw3BMWM3 - 3._dp*TrYb3adjYb3*Yw3adjYw3BMWM3/10._dp - 3._dp*TrYuadjYu*Yw3adjYw3BMWM3 -& 
&  3._dp*TrYw3adjYw3*Yw3adjYw3BMWM3/2._dp - 3._dp*TrYx3adjYx3*Yw3adjYw3BMWM3 -           & 
&  6._dp*g1p2*MassB*Yw3adjYw3MWM3/5._dp + 2._dp*g2p2*MassWB*Yw3adjYw3MWM3 -              & 
&  3._dp*TrTYb3adjYb3*Yw3adjYw3MWM3/5._dp - 6._dp*TrTYuadjYu*Yw3adjYw3MWM3 -             & 
&  3._dp*TrTYw3adjYw3*Yw3adjYw3MWM3 - 6._dp*TrTYx3adjYx3*Yw3adjYw3MWM3 - 3._dp*Yw3adjYw3TYw3adjYw3MWM3 -& 
&  3._dp*Yw3adjYw3Yw3adjYw3BMWM3/2._dp + 40._dp*g2p4*BMWM3 - 64._dp*g2p4*MassWB*MWM3*NGHw3 +& 
&  16._dp*g2p4*BMWM3*NGHw3 - 48._dp*g2p4*MassWB*MWM3*NGHx3  
betaBMWM32 =  betaBMWM32+ 12._dp*g2p4*BMWM3*NGHx3 - 48._dp*g2p4*MassWB*MWM3*NGHxb3 +& 
&  12._dp*g2p4*BMWM3*NGHxb3

 
DBMWM3 = oo16pi2*( betaBMWM31 + oo16pi2 * betaBMWM32 ) 

 
Else 
DBMWM3 = oo16pi2* betaBMWM31 
End If 
 
 
!-------------------- 
! BMGM3 
!-------------------- 
 
betaBMGM31  = 24._dp*g3p2*MassG*MGM3 - 12._dp*g3p2*BMGM3

 
 
If (TwoLoopRGE) Then 
betaBMGM32 = -144._dp*g3p4*MassG*MGM3 + 36._dp*g3p4*BMGM3 - 144._dp*g3p4*MassG*MGM3*NGHg3 +& 
&  36._dp*g3p4*BMGM3*NGHg3 - 48._dp*g3p4*MassG*MGM3*NGHx3 +& 
&  12._dp*g3p4*BMGM3*NGHx3 - 48._dp*g3p4*MassG*MGM3*NGHxb3 +& 
&  12._dp*g3p4*BMGM3*NGHxb3

 
DBMGM3 = oo16pi2*( betaBMGM31 + oo16pi2 * betaBMGM32 ) 

 
Else 
DBMGM3 = oo16pi2* betaBMGM31 
End If 
 
 
!-------------------- 
! BMBM3 
!-------------------- 
 
betaBMBM31  = 3._dp*BMBM3CYb3TpYb3/5._dp + 6._dp*MBM3CYb3TpTYb3/5._dp +               & 
&  6._dp*TYb3adjYb3MBM3/5._dp + 3._dp*Yb3adjYb3BMBM3/5._dp

 
 
If (TwoLoopRGE) Then 
betaBMBM32 = -9._dp*BMBM3CYb3TpYb3CYb3TpYb3/50._dp - 3._dp*BMBM3CYb3TpYeCYeTpYb3/5._dp -           & 
&  9._dp*BMBM3CYb3TpYw3CYw3TpYb3/10._dp + 9._dp*BMBM3CYb3TpYb3*g1p2/25._dp +             & 
&  9._dp*BMBM3CYb3TpYb3*g2p2/5._dp + 18._dp*g1p2*MBM3CYb3TpTYb3/25._dp + 18._dp*g2p2*MBM3CYb3TpTYb3/5._dp -& 
&  9._dp*MBM3CYb3TpTYb3CYb3TpYb3/25._dp - 6._dp*MBM3CYb3TpTYeCYeTpYb3/5._dp -            & 
&  9._dp*MBM3CYb3TpTYw3CYw3TpYb3/5._dp - 18._dp*g1p2*MassB*MBM3CYb3TpYb3/25._dp -        & 
&  18._dp*g2p2*MassWB*MBM3CYb3TpYb3/5._dp - 9._dp*MBM3CYb3TpYb3CYb3TpTYb3/25._dp -       & 
&  6._dp*MBM3CYb3TpYeCYeTpTYb3/5._dp - 9._dp*MBM3CYb3TpYw3CYw3TpTYb3/5._dp -             & 
&  9._dp*MBM3CYb3TpYb3*TrTYb3adjYb3/25._dp - 18._dp*MBM3CYb3TpYb3*TrTYuadjYu/5._dp -     & 
&  9._dp*MBM3CYb3TpYb3*TrTYw3adjYw3/5._dp - 18._dp*MBM3CYb3TpYb3*TrTYx3adjYx3/5._dp -    & 
&  9._dp*BMBM3CYb3TpYb3*TrYb3adjYb3/50._dp - 9._dp*MBM3CYb3TpTYb3*TrYb3adjYb3/25._dp -   & 
&  9._dp*BMBM3CYb3TpYb3*TrYuadjYu/5._dp - 18._dp*MBM3CYb3TpTYb3*TrYuadjYu/5._dp  
betaBMBM32 =  betaBMBM32- 9._dp*BMBM3CYb3TpYb3*TrYw3adjYw3/10._dp - 9._dp*MBM3CYb3TpTYb3*TrYw3adjYw3/5._dp -    & 
&  9._dp*BMBM3CYb3TpYb3*TrYx3adjYx3/5._dp - 18._dp*MBM3CYb3TpTYb3*TrYx3adjYx3/5._dp +    & 
&  18._dp*g1p2*TYb3adjYb3MBM3/25._dp + 18._dp*g2p2*TYb3adjYb3MBM3/5._dp - 9._dp*TrYb3adjYb3*TYb3adjYb3MBM3/25._dp -& 
&  18._dp*TrYuadjYu*TYb3adjYb3MBM3/5._dp - 9._dp*TrYw3adjYw3*TYb3adjYb3MBM3/5._dp -      & 
&  18._dp*TrYx3adjYx3*TYb3adjYb3MBM3/5._dp - 9._dp*TYb3adjYb3Yb3adjYb3MBM3/25._dp -      & 
&  6._dp*TYb3adjYeYeadjYb3MBM3/5._dp - 9._dp*TYb3adjYw3Yw3adjYb3MBM3/5._dp +             & 
&  9._dp*g1p2*Yb3adjYb3BMBM3/25._dp + 9._dp*g2p2*Yb3adjYb3BMBM3/5._dp - 9._dp*TrYb3adjYb3*Yb3adjYb3BMBM3/50._dp -& 
&  9._dp*TrYuadjYu*Yb3adjYb3BMBM3/5._dp - 9._dp*TrYw3adjYw3*Yb3adjYb3BMBM3/10._dp -      & 
&  9._dp*TrYx3adjYx3*Yb3adjYb3BMBM3/5._dp - 18._dp*g1p2*MassB*Yb3adjYb3MBM3/25._dp -     & 
&  18._dp*g2p2*MassWB*Yb3adjYb3MBM3/5._dp - 9._dp*TrTYb3adjYb3*Yb3adjYb3MBM3/25._dp -    & 
&  18._dp*TrTYuadjYu*Yb3adjYb3MBM3/5._dp - 9._dp*TrTYw3adjYw3*Yb3adjYb3MBM3/5._dp  
betaBMBM32 =  betaBMBM32- 18._dp*TrTYx3adjYx3*Yb3adjYb3MBM3/5._dp - 9._dp*Yb3adjYb3TYb3adjYb3MBM3/25._dp -      & 
&  9._dp*Yb3adjYb3Yb3adjYb3BMBM3/50._dp - 6._dp*Yb3adjYeTYeadjYb3MBM3/5._dp -            & 
&  3._dp*Yb3adjYeYeadjYb3BMBM3/5._dp - 9._dp*Yb3adjYw3TYw3adjYb3MBM3/5._dp -             & 
&  9._dp*Yb3adjYw3Yw3adjYb3BMBM3/10._dp

 
DBMBM3 = oo16pi2*( betaBMBM31 + oo16pi2 * betaBMBM32 ) 

 
Else 
DBMBM3 = oo16pi2* betaBMBM31 
End If 
 
 
!-------------------- 
! mq2 
!-------------------- 
 
betamq21  = mq2TpYdCYd + mq2TpYuCYu + 2._dp*TpTYdCTYd + 2._dp*TpTYuCTYu +             & 
&  2._dp*mHd2*TpYdCYd + TpYdCYdmq2 + 2._dp*TpYdmd2CYd + 2._dp*mHu2*TpYuCYu +             & 
&  TpYuCYumq2 + 2._dp*TpYumu2CYu - 2._dp*g1p2*id3R*MassB*Conjg(MassB)/15._dp -           & 
&  32._dp*g3p2*id3R*MassG*Conjg(MassG)/3._dp - 6._dp*g2p2*id3R*MassWB*Conjg(MassWB)      & 
&  + (g1p2*id3R*Tr1(1))/3._dp

 
 
If (TwoLoopRGE) Then 
betamq22 = -4._dp*adjTYdTYdTpYdCYd - 4._dp*adjTYuTYuTpYuCYu + 2._dp*g1p2*mq2TpYdCYd/5._dp -      & 
&  2._dp*mq2TpYdCYdTpYdCYd + 4._dp*g1p2*mq2TpYuCYu/5._dp - 2._dp*mq2TpYuCYuTpYuCYu -     & 
&  4._dp*TpTYdadjYdTpYx3CTYx3 - 4._dp*TpTYdadjYx3Yx3CTYd + 4._dp*g1p2*TpTYdCTYd/5._dp -  & 
&  4._dp*TpTYdCYdTpYdCTYd + 8._dp*g1p2*TpTYuCTYu/5._dp - 4._dp*TpTYuCYuTpYuCTYu -        & 
&  4._dp*TpYdadjYdTpTYx3CTYx3 - 4._dp*TpYdadjYx3mHxb32Yx3CYd - 4._dp*TpYdadjYx3TYx3CTYd -& 
&  4._dp*mHd2*TpYdadjYx3Yx3CYd - 4._dp*mHu2*TpYdadjYx3Yx3CYd - 2._dp*TpYdadjYx3Yx3CYdmq2 -& 
&  4._dp*TpYdadjYx3Yx3md2CYd - 4._dp*g1p2*MassB*TpYdCTYd/5._dp - 4._dp*TpYdCTYdTpTYdCYd +& 
&  4._dp*g1p2*mHd2*TpYdCYd/5._dp + 2._dp*g1p2*TpYdCYdmq2/5._dp - 4._dp*TpYdCYdmq2TpYdCYd -& 
&  4._dp*TpYdCYdTpTYdCTYd - 8._dp*mHd2*TpYdCYdTpYdCYd - 2._dp*TpYdCYdTpYdCYdmq2 -        & 
&  4._dp*TpYdCYdTpYdmd2CYd - 4._dp*TpYdmd2adjYx3Yx3CYd + 4._dp*g1p2*TpYdmd2CYd/5._dp -   & 
&  4._dp*TpYdmd2CYdTpYdCYd - 8._dp*g1p2*MassB*TpYuCTYu/5._dp - 4._dp*TpYuCTYuTpTYuCYu  
betamq22 =  betamq22+ 8._dp*g1p2*mHu2*TpYuCYu/5._dp + 4._dp*g1p2*TpYuCYumq2/5._dp - 4._dp*TpYuCYumq2TpYuCYu -& 
&  4._dp*TpYuCYuTpTYuCTYu - 8._dp*mHu2*TpYuCYuTpYuCYu - 2._dp*TpYuCYuTpYuCYumq2 -        & 
&  4._dp*TpYuCYuTpYumu2CYu + 8._dp*g1p2*TpYumu2CYu/5._dp - 4._dp*TpYumu2CYuTpYuCYu -     & 
&  6._dp*TpYdCYd*Trmd2YdadjYd - 2._dp*TpYdCYd*Trme2YeadjYe - 3._dp*TpYuCYu*TrmHb32Yb3adjYb3/5._dp -& 
&  3._dp*TpYuCYu*TrmHw32Yw3adjYw3 - 6._dp*TpYuCYu*TrmHxb32Yx3adjYx3 - 6._dp*TpYuCYu*Trmu2YuadjYu -& 
&  3._dp*TpYuCYu*TrTpTYb3CTYb3/5._dp - 6._dp*TpYdCYd*TrTpTYdCTYd - 2._dp*TpYdCYd*TrTpTYeCTYe -& 
&  6._dp*TpYuCYu*TrTpTYuCTYu - 3._dp*TpYuCYu*TrTpTYw3CTYw3 - 6._dp*TpYuCYu*TrTpTYx3CTYx3 -& 
&  3._dp*TpTYuCYu*TrTpYb3CTYb3/5._dp - 6._dp*TpTYdCYd*TrTpYdCTYd - 2._dp*TpTYdCYd*TrTpYeCTYe -& 
&  6._dp*TpTYuCYu*TrTpYuCTYu - 3._dp*TpTYuCYu*TrTpYw3CTYw3 - 6._dp*TpTYuCYu*TrTpYx3CTYx3 -& 
&  3._dp*TpYuCTYu*TrTYb3adjYb3/5._dp - 6._dp*TpYdCTYd*TrTYdadjYd - 2._dp*TpYdCTYd*TrTYeadjYe -& 
&  6._dp*TpYuCTYu*TrTYuadjYu - 3._dp*TpYuCTYu*TrTYw3adjYw3 - 6._dp*TpYuCTYu*TrTYx3adjYx3  
betamq22 =  betamq22- 3._dp*mq2TpYuCYu*TrYb3adjYb3/10._dp - 3._dp*TpTYuCTYu*TrYb3adjYb3/5._dp -             & 
&  6._dp*mHu2*TpYuCYu*TrYb3adjYb3/5._dp - 3._dp*TpYuCYumq2*TrYb3adjYb3/10._dp -          & 
&  3._dp*TpYumu2CYu*TrYb3adjYb3/5._dp - 3._dp*TpYuCYu*TrYb3ml2adjYb3/5._dp -             & 
&  3._dp*mq2TpYdCYd*TrYdadjYd - 6._dp*TpTYdCTYd*TrYdadjYd - 12._dp*mHd2*TpYdCYd*TrYdadjYd -& 
&  3._dp*TpYdCYdmq2*TrYdadjYd - 6._dp*TpYdmd2CYd*TrYdadjYd - 6._dp*TpYdCYd*TrYdmq2adjYd -& 
&  mq2TpYdCYd*TrYeadjYe - 2._dp*TpTYdCTYd*TrYeadjYe - 4._dp*mHd2*TpYdCYd*TrYeadjYe -     & 
&  TpYdCYdmq2*TrYeadjYe - 2._dp*TpYdmd2CYd*TrYeadjYe - 2._dp*TpYdCYd*TrYeml2adjYe -      & 
&  3._dp*mq2TpYuCYu*TrYuadjYu - 6._dp*TpTYuCTYu*TrYuadjYu - 12._dp*mHu2*TpYuCYu*TrYuadjYu -& 
&  3._dp*TpYuCYumq2*TrYuadjYu - 6._dp*TpYumu2CYu*TrYuadjYu - 6._dp*TpYuCYu*TrYumq2adjYu -& 
&  3._dp*mq2TpYuCYu*TrYw3adjYw3/2._dp - 3._dp*TpTYuCTYu*TrYw3adjYw3 - 6._dp*mHu2*TpYuCYu*TrYw3adjYw3 -& 
&  3._dp*TpYuCYumq2*TrYw3adjYw3/2._dp - 3._dp*TpYumu2CYu*TrYw3adjYw3 - 3._dp*TpYuCYu*TrYw3ml2adjYw3  
betamq22 =  betamq22- 3._dp*mq2TpYuCYu*TrYx3adjYx3 - 6._dp*TpTYuCTYu*TrYx3adjYx3 - 12._dp*mHu2*TpYuCYu*TrYx3adjYx3 -& 
&  3._dp*TpYuCYumq2*TrYx3adjYx3 - 6._dp*TpYumu2CYu*TrYx3adjYx3 - 6._dp*TpYuCYu*TrYx3md2adjYx3 -& 
&  2._dp*Ydmq2adjYx3Yx3CYd + 199._dp*g1p4*id3R*MassB*Conjg(MassB)/75._dp +               & 
&  2._dp*g1p2*g2p2*id3R*MassB*Conjg(MassB)/5._dp + 32._dp*g1p2*g3p2*id3R*MassB*Conjg(MassB)/45._dp +& 
&  16._dp*g1p2*g3p2*id3R*MassG*Conjg(MassB)/45._dp + (g1p2*g2p2*id3R*MassWB*Conjg(MassB))/5._dp -& 
&  4._dp*g1p2*TpTYdCYd*Conjg(MassB)/5._dp - 8._dp*g1p2*TpTYuCYu*Conjg(MassB)/5._dp +     & 
&  8._dp*g1p2*MassB*TpYdCYd*Conjg(MassB)/5._dp + 16._dp*g1p2*MassB*TpYuCYu*Conjg(MassB)/5._dp +& 
&  16._dp*g1p2*g3p2*id3R*MassB*Conjg(MassG)/45._dp + 32._dp*g1p2*g3p2*id3R*MassG*Conjg(MassG)/45._dp +& 
&  32._dp*g2p2*g3p2*id3R*MassG*Conjg(MassG) - 128._dp*g3p4*id3R*MassG*Conjg(MassG)/3._dp +& 
&  16._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassG) + (g1p2*g2p2*id3R*MassB*Conjg(MassWB))/5._dp +& 
&  16._dp*g2p2*g3p2*id3R*MassG*Conjg(MassWB) + 2._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassWB)/5._dp  
betamq22 =  betamq22+ 33._dp*g2p4*id3R*MassWB*Conjg(MassWB) + 32._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassWB) +  & 
&  96._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 + 36._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHw3 +& 
&  g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 +& 
&  27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHx3 + g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 +& 
&  32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3 + 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHxb3 +& 
&  2._dp*g1p4*id3R*Tr2(1)/15._dp + 6._dp*g2p4*id3R*Tr2(2) + 32._dp*g3p4*id3R*Tr2(3)/3._dp +& 
&  4._dp*g1p2*id3R*Tr3(1)/3._dp

 
Dmq2 = oo16pi2*( betamq21 + oo16pi2 * betamq22 ) 

 
Else 
Dmq2 = oo16pi2* betamq21 
End If 
 
 
Forall(i1=1:3) Dmq2(i1,i1) =  Real(Dmq2(i1,i1),dp) 
!-------------------- 
! ml2 
!-------------------- 
 
betaml21  = 3._dp*ml2TpYb3CYb3/10._dp + ml2TpYeCYe + 3._dp*ml2TpYw3CYw3/2._dp +       & 
&  3._dp*TpTYb3CTYb3/5._dp + 2._dp*TpTYeCTYe + 3._dp*TpTYw3CTYw3 + 3._dp*mHu2*TpYb3CYb3/5._dp +& 
&  3._dp*TpYb3CYb3ml2/10._dp + 3._dp*TpYb3mHb32CYb3/5._dp + 2._dp*mHd2*TpYeCYe +         & 
&  TpYeCYeml2 + 2._dp*TpYeme2CYe + 3._dp*mHu2*TpYw3CYw3 + 3._dp*TpYw3CYw3ml2/2._dp +     & 
&  3._dp*TpYw3mHw32CYw3 - 6._dp*g1p2*id3R*MassB*Conjg(MassB)/5._dp - 6._dp*g2p2*id3R*MassWB*Conjg(MassWB)& 
&  - g1p2*id3R*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betaml22 = -9._dp*adjTYb3TYb3TpYb3CYb3/25._dp - 9._dp*adjTYb3TYb3TpYw3CYw3/20._dp -              & 
&  4._dp*adjTYeTYeTpYeCYe - 9._dp*adjTYw3TYw3TpYb3CYb3/20._dp - 9._dp*adjTYw3TYw3TpYw3CYw3/4._dp -& 
&  9._dp*ml2TpYb3CYb3TpYb3CYb3/50._dp - 9._dp*ml2TpYb3CYb3TpYw3CYw3/40._dp +             & 
&  6._dp*g1p2*ml2TpYeCYe/5._dp - 2._dp*ml2TpYeCYeTpYeCYe + 6._dp*g2p2*ml2TpYw3CYw3 -     & 
&  9._dp*ml2TpYw3CYw3TpYb3CYb3/40._dp - 9._dp*ml2TpYw3CYw3TpYw3CYw3/8._dp -              & 
&  9._dp*TpTYb3CYb3TpYb3CTYb3/25._dp - 9._dp*TpTYb3CYb3TpYw3CTYw3/20._dp +               & 
&  12._dp*g1p2*TpTYeCTYe/5._dp - 4._dp*TpTYeCYeTpYeCTYe + 12._dp*g2p2*TpTYw3CTYw3 -      & 
&  9._dp*TpTYw3CYw3TpYb3CTYb3/20._dp - 9._dp*TpTYw3CYw3TpYw3CTYw3/4._dp - 9._dp*TpYb3CTYb3TpTYb3CYb3/25._dp -& 
&  9._dp*TpYb3CTYb3TpTYw3CYw3/20._dp - 9._dp*TpYb3CYb3ml2TpYb3CYb3/25._dp -              & 
&  9._dp*TpYb3CYb3ml2TpYw3CYw3/20._dp - 9._dp*TpYb3CYb3TpTYb3CTYb3/25._dp -              & 
&  9._dp*TpYb3CYb3TpTYw3CTYw3/20._dp - 18._dp*mHu2*TpYb3CYb3TpYb3CYb3/25._dp  
betaml22 =  betaml22- 9._dp*TpYb3CYb3TpYb3CYb3ml2/50._dp - 9._dp*TpYb3CYb3TpYb3mHb32CYb3/25._dp -           & 
&  9._dp*mHu2*TpYb3CYb3TpYw3CYw3/10._dp - 9._dp*TpYb3CYb3TpYw3CYw3ml2/40._dp -           & 
&  9._dp*TpYb3CYb3TpYw3mHw32CYw3/20._dp - 9._dp*TpYb3mHb32CYb3TpYb3CYb3/25._dp -         & 
&  9._dp*TpYb3mHb32CYb3TpYw3CYw3/20._dp - 12._dp*g1p2*MassB*TpYeCTYe/5._dp -             & 
&  4._dp*TpYeCTYeTpTYeCYe + 12._dp*g1p2*mHd2*TpYeCYe/5._dp + 6._dp*g1p2*TpYeCYeml2/5._dp -& 
&  4._dp*TpYeCYeml2TpYeCYe - 4._dp*TpYeCYeTpTYeCTYe - 8._dp*mHd2*TpYeCYeTpYeCYe -        & 
&  2._dp*TpYeCYeTpYeCYeml2 - 4._dp*TpYeCYeTpYeme2CYe + 12._dp*g1p2*TpYeme2CYe/5._dp -    & 
&  4._dp*TpYeme2CYeTpYeCYe - 12._dp*g2p2*MassWB*TpYw3CTYw3 - 9._dp*TpYw3CTYw3TpTYb3CYb3/20._dp -& 
&  9._dp*TpYw3CTYw3TpTYw3CYw3/4._dp + 12._dp*g2p2*mHu2*TpYw3CYw3 + 6._dp*g2p2*TpYw3CYw3ml2 -& 
&  9._dp*TpYw3CYw3ml2TpYb3CYb3/20._dp - 9._dp*TpYw3CYw3ml2TpYw3CYw3/4._dp -              & 
&  9._dp*TpYw3CYw3TpTYb3CTYb3/20._dp - 9._dp*TpYw3CYw3TpTYw3CTYw3/4._dp - 9._dp*mHu2*TpYw3CYw3TpYb3CYb3/10._dp  
betaml22 =  betaml22- 9._dp*TpYw3CYw3TpYb3CYb3ml2/40._dp - 9._dp*TpYw3CYw3TpYb3mHb32CYb3/20._dp -           & 
&  9._dp*mHu2*TpYw3CYw3TpYw3CYw3/2._dp - 9._dp*TpYw3CYw3TpYw3CYw3ml2/8._dp -             & 
&  9._dp*TpYw3CYw3TpYw3mHw32CYw3/4._dp + 12._dp*g2p2*TpYw3mHw32CYw3 - 9._dp*TpYw3mHw32CYw3TpYb3CYb3/20._dp -& 
&  9._dp*TpYw3mHw32CYw3TpYw3CYw3/4._dp - 6._dp*TpYeCYe*Trmd2YdadjYd - 2._dp*TpYeCYe*Trme2YeadjYe -& 
&  9._dp*TpYb3CYb3*TrmHb32Yb3adjYb3/50._dp - 9._dp*TpYw3CYw3*TrmHb32Yb3adjYb3/10._dp -   & 
&  9._dp*TpYb3CYb3*TrmHw32Yw3adjYw3/10._dp - 9._dp*TpYw3CYw3*TrmHw32Yw3adjYw3/2._dp -    & 
&  9._dp*TpYb3CYb3*TrmHxb32Yx3adjYx3/5._dp - 9._dp*TpYw3CYw3*TrmHxb32Yx3adjYx3 -         & 
&  9._dp*TpYb3CYb3*Trmu2YuadjYu/5._dp - 9._dp*TpYw3CYw3*Trmu2YuadjYu - 9._dp*TpYb3CYb3*TrTpTYb3CTYb3/50._dp -& 
&  9._dp*TpYw3CYw3*TrTpTYb3CTYb3/10._dp - 6._dp*TpYeCYe*TrTpTYdCTYd - 2._dp*TpYeCYe*TrTpTYeCTYe -& 
&  9._dp*TpYb3CYb3*TrTpTYuCTYu/5._dp - 9._dp*TpYw3CYw3*TrTpTYuCTYu - 9._dp*TpYb3CYb3*TrTpTYw3CTYw3/10._dp -& 
&  9._dp*TpYw3CYw3*TrTpTYw3CTYw3/2._dp - 9._dp*TpYb3CYb3*TrTpTYx3CTYx3/5._dp  
betaml22 =  betaml22- 9._dp*TpYw3CYw3*TrTpTYx3CTYx3 - 9._dp*TpTYb3CYb3*TrTpYb3CTYb3/50._dp - 9._dp*TpTYw3CYw3*TrTpYb3CTYb3/10._dp -& 
&  6._dp*TpTYeCYe*TrTpYdCTYd - 2._dp*TpTYeCYe*TrTpYeCTYe - 9._dp*TpTYb3CYb3*TrTpYuCTYu/5._dp -& 
&  9._dp*TpTYw3CYw3*TrTpYuCTYu - 9._dp*TpTYb3CYb3*TrTpYw3CTYw3/10._dp - 9._dp*TpTYw3CYw3*TrTpYw3CTYw3/2._dp -& 
&  9._dp*TpTYb3CYb3*TrTpYx3CTYx3/5._dp - 9._dp*TpTYw3CYw3*TrTpYx3CTYx3 - 9._dp*TpYb3CTYb3*TrTYb3adjYb3/50._dp -& 
&  9._dp*TpYw3CTYw3*TrTYb3adjYb3/10._dp - 6._dp*TpYeCTYe*TrTYdadjYd - 2._dp*TpYeCTYe*TrTYeadjYe -& 
&  9._dp*TpYb3CTYb3*TrTYuadjYu/5._dp - 9._dp*TpYw3CTYw3*TrTYuadjYu - 9._dp*TpYb3CTYb3*TrTYw3adjYw3/10._dp -& 
&  9._dp*TpYw3CTYw3*TrTYw3adjYw3/2._dp - 9._dp*TpYb3CTYb3*TrTYx3adjYx3/5._dp -           & 
&  9._dp*TpYw3CTYw3*TrTYx3adjYx3 - 9._dp*ml2TpYb3CYb3*TrYb3adjYb3/100._dp -              & 
&  9._dp*ml2TpYw3CYw3*TrYb3adjYb3/20._dp - 9._dp*TpTYb3CTYb3*TrYb3adjYb3/50._dp -        & 
&  9._dp*TpTYw3CTYw3*TrYb3adjYb3/10._dp - 9._dp*mHu2*TpYb3CYb3*TrYb3adjYb3/25._dp -      & 
&  9._dp*TpYb3CYb3ml2*TrYb3adjYb3/100._dp - 9._dp*TpYb3mHb32CYb3*TrYb3adjYb3/50._dp  
betaml22 =  betaml22- 9._dp*mHu2*TpYw3CYw3*TrYb3adjYb3/5._dp - 9._dp*TpYw3CYw3ml2*TrYb3adjYb3/20._dp -      & 
&  9._dp*TpYw3mHw32CYw3*TrYb3adjYb3/10._dp - 9._dp*TpYb3CYb3*TrYb3ml2adjYb3/50._dp -     & 
&  9._dp*TpYw3CYw3*TrYb3ml2adjYb3/10._dp - 3._dp*ml2TpYeCYe*TrYdadjYd - 6._dp*TpTYeCTYe*TrYdadjYd -& 
&  12._dp*mHd2*TpYeCYe*TrYdadjYd - 3._dp*TpYeCYeml2*TrYdadjYd - 6._dp*TpYeme2CYe*TrYdadjYd -& 
&  6._dp*TpYeCYe*TrYdmq2adjYd - ml2TpYeCYe*TrYeadjYe - 2._dp*TpTYeCTYe*TrYeadjYe -       & 
&  4._dp*mHd2*TpYeCYe*TrYeadjYe - TpYeCYeml2*TrYeadjYe - 2._dp*TpYeme2CYe*TrYeadjYe -    & 
&  2._dp*TpYeCYe*TrYeml2adjYe - 9._dp*ml2TpYb3CYb3*TrYuadjYu/10._dp - 9._dp*ml2TpYw3CYw3*TrYuadjYu/2._dp -& 
&  9._dp*TpTYb3CTYb3*TrYuadjYu/5._dp - 9._dp*TpTYw3CTYw3*TrYuadjYu - 18._dp*mHu2*TpYb3CYb3*TrYuadjYu/5._dp -& 
&  9._dp*TpYb3CYb3ml2*TrYuadjYu/10._dp - 9._dp*TpYb3mHb32CYb3*TrYuadjYu/5._dp -          & 
&  18._dp*mHu2*TpYw3CYw3*TrYuadjYu - 9._dp*TpYw3CYw3ml2*TrYuadjYu/2._dp - 9._dp*TpYw3mHw32CYw3*TrYuadjYu -& 
&  9._dp*TpYb3CYb3*TrYumq2adjYu/5._dp - 9._dp*TpYw3CYw3*TrYumq2adjYu - 9._dp*ml2TpYb3CYb3*TrYw3adjYw3/20._dp  
betaml22 =  betaml22- 9._dp*ml2TpYw3CYw3*TrYw3adjYw3/4._dp - 9._dp*TpTYb3CTYb3*TrYw3adjYw3/10._dp -         & 
&  9._dp*TpTYw3CTYw3*TrYw3adjYw3/2._dp - 9._dp*mHu2*TpYb3CYb3*TrYw3adjYw3/5._dp -        & 
&  9._dp*TpYb3CYb3ml2*TrYw3adjYw3/20._dp - 9._dp*TpYb3mHb32CYb3*TrYw3adjYw3/10._dp -     & 
&  9._dp*mHu2*TpYw3CYw3*TrYw3adjYw3 - 9._dp*TpYw3CYw3ml2*TrYw3adjYw3/4._dp -             & 
&  9._dp*TpYw3mHw32CYw3*TrYw3adjYw3/2._dp - 9._dp*TpYb3CYb3*TrYw3ml2adjYw3/10._dp -      & 
&  9._dp*TpYw3CYw3*TrYw3ml2adjYw3/2._dp - 9._dp*ml2TpYb3CYb3*TrYx3adjYx3/10._dp -        & 
&  9._dp*ml2TpYw3CYw3*TrYx3adjYx3/2._dp - 9._dp*TpTYb3CTYb3*TrYx3adjYx3/5._dp -          & 
&  9._dp*TpTYw3CTYw3*TrYx3adjYx3 - 18._dp*mHu2*TpYb3CYb3*TrYx3adjYx3/5._dp -             & 
&  9._dp*TpYb3CYb3ml2*TrYx3adjYx3/10._dp - 9._dp*TpYb3mHb32CYb3*TrYx3adjYx3/5._dp -      & 
&  18._dp*mHu2*TpYw3CYw3*TrYx3adjYx3 - 9._dp*TpYw3CYw3ml2*TrYx3adjYx3/2._dp -            & 
&  9._dp*TpYw3mHw32CYw3*TrYx3adjYx3 - 9._dp*TpYb3CYb3*TrYx3md2adjYx3/5._dp  
betaml22 =  betaml22- 9._dp*TpYw3CYw3*TrYx3md2adjYx3 + 621._dp*g1p4*id3R*MassB*Conjg(MassB)/25._dp +        & 
&  18._dp*g1p2*g2p2*id3R*MassB*Conjg(MassB)/5._dp + 9._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassB)/5._dp -& 
&  12._dp*g1p2*TpTYeCYe*Conjg(MassB)/5._dp + 24._dp*g1p2*MassB*TpYeCYe*Conjg(MassB)/5._dp +& 
&  9._dp*g1p2*g2p2*id3R*MassB*Conjg(MassWB)/5._dp + 18._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassWB)/5._dp +& 
&  33._dp*g2p4*id3R*MassWB*Conjg(MassWB) - 12._dp*g2p2*TpTYw3CYw3*Conjg(MassWB) +        & 
&  24._dp*g2p2*MassWB*TpYw3CYw3*Conjg(MassWB) + 36._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHw3 +& 
&  9._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHx3 +& 
&  9._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 + 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHxb3 +& 
&  6._dp*g1p4*id3R*Tr2(1)/5._dp + 6._dp*g2p4*id3R*Tr2(2) - 4._dp*g1p2*id3R*Tr3(1)

 
Dml2 = oo16pi2*( betaml21 + oo16pi2 * betaml22 ) 

 
Else 
Dml2 = oo16pi2* betaml21 
End If 
 
 
Forall(i1=1:3) Dml2(i1,i1) =  Real(Dml2(i1,i1),dp) 
!-------------------- 
! mHd2 
!-------------------- 
 
betamHd21  = 6._dp*Trmd2YdadjYd + 2._dp*Trme2YeadjYe + 6._dp*TrTpTYdCTYd +            & 
&  2._dp*TrTpTYeCTYe + 6._dp*mHd2*TrYdadjYd + 6._dp*TrYdmq2adjYd + 2._dp*mHd2*TrYeadjYe +& 
&  2._dp*TrYeml2adjYe - 6._dp*g1p2*MassB*Conjg(MassB)/5._dp - 6._dp*g2p2*MassWB*Conjg(MassWB)& 
&  - g1p2*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHd22 = -12._dp*TrCTYdTpYdadjYx3TYx3 - 12._dp*TrCYdTpYdadjYx3mHxb32Yx3 - 4._dp*g1p2*Trmd2YdadjYd/5._dp +& 
&  32._dp*g3p2*Trmd2YdadjYd - 12._dp*Trmd2YdadjYdTpYx3CYx3 + 12._dp*g1p2*Trme2YeadjYe/5._dp -& 
&  12._dp*TrTpTYdadjYx3Yx3CTYd - 4._dp*g1p2*TrTpTYdCTYd/5._dp + 32._dp*g3p2*TrTpTYdCTYd +& 
&  12._dp*g1p2*TrTpTYeCTYe/5._dp + 4._dp*g1p2*MassB*TrTpYdCTYd/5._dp - 32._dp*g3p2*MassG*TrTpYdCTYd -& 
&  12._dp*g1p2*MassB*TrTpYeCTYe/5._dp - 3._dp*TrTYb3adjYeYeadjTYb3/5._dp -               & 
&  12._dp*TrTYdadjYdTpYx3CTYx3 - 18._dp*TrTYdadjYdYdadjTYd - 6._dp*TrTYeadjYeYeadjTYe -  & 
&  6._dp*TrTYuadjYdYdadjTYu - 3._dp*TrTYw3adjYeYeadjTYw3 - 3._dp*TrYb3adjYeme2YeadjYb3/5._dp -& 
&  3._dp*TrYb3adjYeTYeadjTYb3/5._dp - 3._dp*mHd2*TrYb3adjYeYeadjYb3/10._dp -             & 
&  3._dp*mHu2*TrYb3adjYeYeadjYb3/5._dp - 3._dp*TrYb3adjYeYeml2adjYb3/5._dp -             & 
&  3._dp*TrYb3ml2adjYeYeadjYb3/5._dp - 3._dp*TrYb3TpTYeCTYeadjYb3/5._dp - 4._dp*g1p2*mHd2*TrYdadjYd/5._dp +& 
&  32._dp*g3p2*mHd2*TrYdadjYd - 12._dp*TrYdadjYdmd2TpYx3CYx3 - 36._dp*TrYdadjYdmd2YdadjYd  
betamHd22 =  betamHd22- 12._dp*TrYdadjYdTpTYx3CTYx3 - 12._dp*mHd2*TrYdadjYdTpYx3CYx3 - 12._dp*mHu2*TrYdadjYdTpYx3CYx3 -& 
&  36._dp*TrYdadjYdTYdadjTYd - 36._dp*mHd2*TrYdadjYdYdadjYd - 18._dp*TrYdadjYdYdmq2adjYd -& 
&  6._dp*TrYdadjYumu2YuadjYd - 6._dp*TrYdadjYuTYuadjTYd - 3._dp*mHd2*TrYdadjYuYuadjYd -  & 
&  4._dp*g1p2*TrYdmq2adjYd/5._dp + 32._dp*g3p2*TrYdmq2adjYd - 12._dp*TrYdmq2adjYdTpYx3CYx3 -& 
&  18._dp*TrYdmq2adjYdYdadjYd - 18._dp*TrYdTpTYdCTYdadjYd - 3._dp*TrYeadjYb3mHb32Yb3adjYe/5._dp -& 
&  3._dp*TrYeadjYb3TYb3adjTYe/5._dp - 3._dp*mHd2*TrYeadjYb3Yb3adjYe/10._dp +             & 
&  12._dp*g1p2*mHd2*TrYeadjYe/5._dp - 12._dp*TrYeadjYeme2YeadjYe - 12._dp*TrYeadjYeTYeadjTYe -& 
&  12._dp*mHd2*TrYeadjYeYeadjYe - 6._dp*TrYeadjYeYeml2adjYe - 3._dp*TrYeadjYw3mHw32Yw3adjYe -& 
&  3._dp*TrYeadjYw3TYw3adjTYe - 3._dp*mHd2*TrYeadjYw3Yw3adjYe/2._dp + 12._dp*g1p2*TrYeml2adjYe/5._dp -& 
&  6._dp*TrYeml2adjYeYeadjYe - 6._dp*TrYeTpTYeCTYeadjYe - 6._dp*TrYuadjYdmd2YdadjYu -    & 
&  6._dp*TrYuadjYdTYdadjTYu - 3._dp*mHd2*TrYuadjYdYdadjYu - 6._dp*mHu2*TrYuadjYdYdadjYu  
betamHd22 =  betamHd22- 6._dp*TrYuadjYdYdmq2adjYu - 6._dp*TrYumq2adjYdYdadjYu - 6._dp*TrYuTpTYdCTYdadjYu -    & 
&  3._dp*TrYw3adjYeme2YeadjYw3 - 3._dp*TrYw3adjYeTYeadjTYw3 - 3._dp*mHd2*TrYw3adjYeYeadjYw3/2._dp -& 
&  3._dp*mHu2*TrYw3adjYeYeadjYw3 - 3._dp*TrYw3adjYeYeml2adjYw3 - 3._dp*TrYw3ml2adjYeYeadjYw3 -& 
&  3._dp*TrYw3TpTYeCTYeadjYw3 + 621._dp*g1p4*MassB*Conjg(MassB)/25._dp + 18._dp*g1p2*g2p2*MassB*Conjg(MassB)/5._dp +& 
&  9._dp*g1p2*g2p2*MassWB*Conjg(MassB)/5._dp + 4._dp*g1p2*TrTYdadjYd*Conjg(MassB)/5._dp -& 
&  12._dp*g1p2*TrTYeadjYe*Conjg(MassB)/5._dp - 8._dp*g1p2*MassB*TrYdadjYd*Conjg(MassB)/5._dp +& 
&  24._dp*g1p2*MassB*TrYeadjYe*Conjg(MassB)/5._dp - 32._dp*g3p2*TrTYdadjYd*Conjg(MassG) +& 
&  64._dp*g3p2*MassG*TrYdadjYd*Conjg(MassG) + 9._dp*g1p2*g2p2*MassB*Conjg(MassWB)/5._dp +& 
&  18._dp*g1p2*g2p2*MassWB*Conjg(MassWB)/5._dp + 33._dp*g2p4*MassWB*Conjg(MassWB) +      & 
&  36._dp*g2p4*MassWB*Conjg(MassWB)*NGHw3 + 9._dp*g1p4*MassB*Conjg(MassB)*NGHx3 +& 
&  27._dp*g2p4*MassWB*Conjg(MassWB)*NGHx3 + 9._dp*g1p4*MassB*Conjg(MassB)*NGHxb3  
betamHd22 =  betamHd22+ 27._dp*g2p4*MassWB*Conjg(MassWB)*NGHxb3 + 6._dp*g1p4*Tr2(1)/5._dp +    & 
&  6._dp*g2p4*Tr2(2) - 4._dp*g1p2*Tr3(1)

 
DmHd2 = oo16pi2*( betamHd21 + oo16pi2 * betamHd22 ) 

 
Else 
DmHd2 = oo16pi2* betamHd21 
End If 
 
 
!-------------------- 
! mHu2 
!-------------------- 
 
betamHu21  = 3._dp*TrmHb32Yb3adjYb3/5._dp + 3._dp*TrmHw32Yw3adjYw3 + 6._dp*TrmHxb32Yx3adjYx3 +& 
&  6._dp*Trmu2YuadjYu + 3._dp*TrTpTYb3CTYb3/5._dp + 6._dp*TrTpTYuCTYu + 3._dp*TrTpTYw3CTYw3 +& 
&  6._dp*TrTpTYx3CTYx3 + 3._dp*mHu2*TrYb3adjYb3/5._dp + 3._dp*TrYb3ml2adjYb3/5._dp +     & 
&  6._dp*mHu2*TrYuadjYu + 6._dp*TrYumq2adjYu + 3._dp*mHu2*TrYw3adjYw3 + 3._dp*TrYw3ml2adjYw3 +& 
&  6._dp*mHu2*TrYx3adjYx3 + 6._dp*TrYx3md2adjYx3 - 6._dp*g1p2*MassB*Conjg(MassB)         & 
& /5._dp - 6._dp*g2p2*MassWB*Conjg(MassWB) + g1p2*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betamHu22 = -12._dp*Trmd2TpYx3CYx3YdadjYd + 12._dp*g2p2*TrmHw32Yw3adjYw3 + 4._dp*g1p2*TrmHxb32Yx3adjYx3 +& 
&  32._dp*g3p2*TrmHxb32Yx3adjYx3 - 12._dp*Trmq2adjYdTpYx3CYx3Yd + 8._dp*g1p2*Trmu2YuadjYu/5._dp +& 
&  32._dp*g3p2*Trmu2YuadjYu - 12._dp*TrTpTYdadjYdTpYx3CTYx3 - 12._dp*TrTpTYdadjYx3Yx3CTYd +& 
&  8._dp*g1p2*TrTpTYuCTYu/5._dp + 32._dp*g3p2*TrTpTYuCTYu + 12._dp*g2p2*TrTpTYw3CTYw3 +  & 
&  4._dp*g1p2*TrTpTYx3CTYx3 + 32._dp*g3p2*TrTpTYx3CTYx3 - 12._dp*TrTpYdadjYdTpTYx3CTYx3 -& 
&  12._dp*TrTpYdadjYx3TYx3CTYd - 8._dp*g1p2*MassB*TrTpYuCTYu/5._dp - 32._dp*g3p2*MassG*TrTpYuCTYu -& 
&  12._dp*g2p2*MassWB*TrTpYw3CTYw3 - 4._dp*g1p2*MassB*TrTpYx3CTYx3 - 32._dp*g3p2*MassG*TrTpYx3CTYx3 -& 
&  12._dp*TrTpYx3CYx3md2YdadjYd - 12._dp*mHd2*TrTpYx3CYx3YdadjYd - 6._dp*mHu2*TrTpYx3CYx3YdadjYd -& 
&  12._dp*TrTpYx3mHxb32CYx3YdadjYd - 27._dp*TrTYb3adjYb3Yb3adjTYb3/50._dp -              & 
&  27._dp*TrTYb3adjYw3Yw3adjTYb3/20._dp - 6._dp*TrTYdadjYuYuadjTYd - 3._dp*TrTYeadjYb3Yb3adjTYe/5._dp -& 
&  3._dp*TrTYeadjYw3Yw3adjTYe - 18._dp*TrTYuadjYuYuadjTYu - 27._dp*TrTYw3adjYb3Yb3adjTYw3/20._dp  
betamHu22 =  betamHu22-27._dp*TrTYw3adjYw3Yw3adjTYw3/4._dp -18._dp*TrTYx3adjYx3Yx3adjTYx3-27._dp*TrYb3adjYb3mHb32Yb3adjYb3/25._dp-&
&  27._dp*TrYb3adjYb3TYb3adjTYb3/25._dp - 27._dp*mHu2*TrYb3adjYb3Yb3adjYb3/25._dp -      & 
&  27._dp*TrYb3adjYb3Yb3ml2adjYb3/50._dp - 3._dp*TrYb3adjYeme2YeadjYb3/5._dp -           & 
&  3._dp*TrYb3adjYeTYeadjTYb3/5._dp - 3._dp*mHu2*TrYb3adjYeYeadjYb3/10._dp -             & 
&  27._dp*TrYb3adjYw3mHw32Yw3adjYb3/10._dp - 27._dp*TrYb3adjYw3TYw3adjTYb3/10._dp -      & 
&  27._dp*mHu2*TrYb3adjYw3Yw3adjYb3/10._dp - 27._dp*TrYb3adjYw3Yw3ml2adjYb3/20._dp -     & 
&  27._dp*TrYb3ml2adjYb3Yb3adjYb3/50._dp - 27._dp*TrYb3ml2adjYw3Yw3adjYb3/20._dp -       & 
&  27._dp*TrYb3TpTYb3CTYb3adjYb3/50._dp - 27._dp*TrYb3TpTYw3CTYw3adjYb3/20._dp -         & 
&  6._dp*mHu2*TrYdadjYdTpYx3CYx3 - 6._dp*TrYdadjYumu2YuadjYd - 6._dp*TrYdadjYuTYuadjTYd -& 
&  6._dp*mHd2*TrYdadjYuYuadjYd - 3._dp*mHu2*TrYdadjYuYuadjYd - 6._dp*TrYdadjYuYumq2adjYd -& 
&  6._dp*TrYdmq2adjYuYuadjYd - 6._dp*TrYdTpTYuCTYuadjYd - 3._dp*TrYeadjYb3mHb32Yb3adjYe/5._dp  
betamHu22 =  betamHu22- 3._dp*TrYeadjYb3TYb3adjTYe/5._dp - 3._dp*mHd2*TrYeadjYb3Yb3adjYe/5._dp -              & 
&  3._dp*mHu2*TrYeadjYb3Yb3adjYe/10._dp - 3._dp*TrYeadjYb3Yb3ml2adjYe/5._dp -            & 
&  3._dp*TrYeadjYw3mHw32Yw3adjYe - 3._dp*TrYeadjYw3TYw3adjTYe - 3._dp*mHd2*TrYeadjYw3Yw3adjYe -& 
&  3._dp*mHu2*TrYeadjYw3Yw3adjYe/2._dp - 3._dp*TrYeadjYw3Yw3ml2adjYe - 3._dp*TrYeml2adjYb3Yb3adjYe/5._dp -& 
&  3._dp*TrYeml2adjYw3Yw3adjYe - 3._dp*TrYeTpTYb3CTYb3adjYe/5._dp - 3._dp*TrYeTpTYw3CTYw3adjYe -& 
&  6._dp*TrYuadjYdmd2YdadjYu - 6._dp*TrYuadjYdTYdadjTYu - 3._dp*mHu2*TrYuadjYdYdadjYu +  & 
&  8._dp*g1p2*mHu2*TrYuadjYu/5._dp + 32._dp*g3p2*mHu2*TrYuadjYu - 36._dp*TrYuadjYumu2YuadjYu -& 
&  36._dp*TrYuadjYuTYuadjTYu - 36._dp*mHu2*TrYuadjYuYuadjYu - 18._dp*TrYuadjYuYumq2adjYu +& 
&  8._dp*g1p2*TrYumq2adjYu/5._dp + 32._dp*g3p2*TrYumq2adjYu - 18._dp*TrYumq2adjYuYuadjYu -& 
&  18._dp*TrYuTpTYuCTYuadjYu - 27._dp*TrYw3adjYb3mHb32Yb3adjYw3/10._dp - 27._dp*TrYw3adjYb3TYb3adjTYw3/10._dp -& 
&  27._dp*mHu2*TrYw3adjYb3Yb3adjYw3/10._dp - 27._dp*TrYw3adjYb3Yb3ml2adjYw3/20._dp  
betamHu22 =  betamHu22- 3._dp*TrYw3adjYeme2YeadjYw3 - 3._dp*TrYw3adjYeTYeadjTYw3 - 3._dp*mHu2*TrYw3adjYeYeadjYw3/2._dp +& 
&  12._dp*g2p2*mHu2*TrYw3adjYw3 - 27._dp*TrYw3adjYw3mHw32Yw3adjYw3/2._dp -               & 
&  27._dp*TrYw3adjYw3TYw3adjTYw3/2._dp - 27._dp*mHu2*TrYw3adjYw3Yw3adjYw3/2._dp -        & 
&  27._dp*TrYw3adjYw3Yw3ml2adjYw3/4._dp - 27._dp*TrYw3ml2adjYb3Yb3adjYw3/20._dp +        & 
&  12._dp*g2p2*TrYw3ml2adjYw3 - 27._dp*TrYw3ml2adjYw3Yw3adjYw3/4._dp - 27._dp*TrYw3TpTYb3CTYb3adjYw3/20._dp -& 
&  27._dp*TrYw3TpTYw3CTYw3adjYw3/4._dp + 4._dp*g1p2*mHu2*TrYx3adjYx3 + 32._dp*g3p2*mHu2*TrYx3adjYx3 -& 
&  36._dp*TrYx3adjYx3mHxb32Yx3adjYx3 - 36._dp*TrYx3adjYx3TYx3adjTYx3 - 36._dp*mHu2*TrYx3adjYx3Yx3adjYx3 -& 
&  18._dp*TrYx3adjYx3Yx3md2adjYx3 + 4._dp*g1p2*TrYx3md2adjYx3 + 32._dp*g3p2*TrYx3md2adjYx3 -& 
&  18._dp*TrYx3md2adjYx3Yx3adjYx3 - 18._dp*TrYx3TpTYx3CTYx3adjYx3 + 621._dp*g1p4*MassB*Conjg(MassB)/25._dp +& 
&  18._dp*g1p2*g2p2*MassB*Conjg(MassB)/5._dp + 9._dp*g1p2*g2p2*MassWB*Conjg(MassB)/5._dp -& 
&  8._dp*g1p2*TrTYuadjYu*Conjg(MassB)/5._dp - 4._dp*g1p2*TrTYx3adjYx3*Conjg(MassB)  
betamHu22 =  betamHu22+ 16._dp*g1p2*MassB*TrYuadjYu*Conjg(MassB)/5._dp + 8._dp*g1p2*MassB*TrYx3adjYx3*Conjg(MassB) -& 
&  32._dp*g3p2*TrTYuadjYu*Conjg(MassG) - 32._dp*g3p2*TrTYx3adjYx3*Conjg(MassG) +         & 
&  64._dp*g3p2*MassG*TrYuadjYu*Conjg(MassG) + 64._dp*g3p2*MassG*TrYx3adjYx3*Conjg(MassG) +& 
&  9._dp*g1p2*g2p2*MassB*Conjg(MassWB)/5._dp + 18._dp*g1p2*g2p2*MassWB*Conjg(MassWB)/5._dp +& 
&  33._dp*g2p4*MassWB*Conjg(MassWB) - 12._dp*g2p2*TrTYw3adjYw3*Conjg(MassWB) +           & 
&  24._dp*g2p2*MassWB*TrYw3adjYw3*Conjg(MassWB) + 36._dp*g2p4*MassWB*Conjg(MassWB)*NGHw3 +& 
&  9._dp*g1p4*MassB*Conjg(MassB)*NGHx3 + 27._dp*g2p4*MassWB*Conjg(MassWB)*NGHx3 +& 
&  9._dp*g1p4*MassB*Conjg(MassB)*NGHxb3 + 27._dp*g2p4*MassWB*Conjg(MassWB)*NGHxb3 +& 
&  6._dp*g1p4*Tr2(1)/5._dp + 6._dp*g2p4*Tr2(2) + 4._dp*g1p2*Tr3(1)

 
DmHu2 = oo16pi2*( betamHu21 + oo16pi2 * betamHu22 ) 

 
Else 
DmHu2 = oo16pi2* betamHu21 
End If 
 
 
!-------------------- 
! md2 
!-------------------- 
 
betamd21  = 2._dp*md2TpYx3CYx3 + 2._dp*md2YdadjYd + 4._dp*TpTYx3CTYx3 +               & 
&  4._dp*mHu2*TpYx3CYx3 + 2._dp*TpYx3CYx3md2 + 4._dp*TpYx3mHxb32CYx3 + 4._dp*TYdadjTYd + & 
&  4._dp*mHd2*YdadjYd + 2._dp*YdadjYdmd2 + 4._dp*Ydmq2adjYd - 8._dp*g1p2*id3R*MassB*Conjg(MassB)& 
& /15._dp - 32._dp*g3p2*id3R*MassG*Conjg(MassG)/3._dp + 2._dp*g1p2*id3R*Tr1(1)/3._dp

 
 
If (TwoLoopRGE) Then 
betamd22 = -4._dp*adjTYx3TYx3TpYx3CYx3 + 2._dp*g1p2*md2TpYx3CYx3 + 6._dp*g2p2*md2TpYx3CYx3 -     & 
&  2._dp*md2TpYx3CYx3TpYx3CYx3 + 2._dp*g1p2*md2YdadjYd/5._dp + 6._dp*g2p2*md2YdadjYd -   & 
&  2._dp*md2YdadjYdYdadjYd - 2._dp*md2YdadjYuYuadjYd + 4._dp*g1p2*TpTYx3CTYx3 +          & 
&  12._dp*g2p2*TpTYx3CTYx3 - 4._dp*TpTYx3CYx3TpYx3CTYx3 - 4._dp*g1p2*MassB*TpYx3CTYx3 -  & 
&  12._dp*g2p2*MassWB*TpYx3CTYx3 - 4._dp*TpYx3CTYx3TpTYx3CYx3 + 4._dp*g1p2*mHu2*TpYx3CYx3 +& 
&  12._dp*g2p2*mHu2*TpYx3CYx3 + 2._dp*g1p2*TpYx3CYx3md2 + 6._dp*g2p2*TpYx3CYx3md2 -      & 
&  4._dp*TpYx3CYx3md2TpYx3CYx3 - 4._dp*TpYx3CYx3TpTYx3CTYx3 - 8._dp*mHu2*TpYx3CYx3TpYx3CYx3 -& 
&  2._dp*TpYx3CYx3TpYx3CYx3md2 - 4._dp*TpYx3CYx3TpYx3mHxb32CYx3 + 4._dp*g1p2*TpYx3mHxb32CYx3 +& 
&  12._dp*g2p2*TpYx3mHxb32CYx3 - 4._dp*TpYx3mHxb32CYx3TpYx3CYx3 - 6._dp*TpYx3CYx3*TrmHb32Yb3adjYb3/5._dp -& 
&  6._dp*TpYx3CYx3*TrmHw32Yw3adjYw3 - 12._dp*TpYx3CYx3*TrmHxb32Yx3adjYx3 -               & 
&  12._dp*TpYx3CYx3*Trmu2YuadjYu - 6._dp*TpYx3CYx3*TrTpTYb3CTYb3/5._dp - 12._dp*TpYx3CYx3*TrTpTYuCTYu  
betamd22 =  betamd22- 6._dp*TpYx3CYx3*TrTpTYw3CTYw3 - 12._dp*TpYx3CYx3*TrTpTYx3CTYx3 - 6._dp*TpTYx3CYx3*TrTpYb3CTYb3/5._dp -& 
&  12._dp*TpTYx3CYx3*TrTpYuCTYu - 6._dp*TpTYx3CYx3*TrTpYw3CTYw3 - 12._dp*TpTYx3CYx3*TrTpYx3CTYx3 -& 
&  6._dp*TpYx3CTYx3*TrTYb3adjYb3/5._dp - 12._dp*TpYx3CTYx3*TrTYuadjYu - 6._dp*TpYx3CTYx3*TrTYw3adjYw3 -& 
&  12._dp*TpYx3CTYx3*TrTYx3adjYx3 - 3._dp*md2TpYx3CYx3*TrYb3adjYb3/5._dp -               & 
&  6._dp*TpTYx3CTYx3*TrYb3adjYb3/5._dp - 12._dp*mHu2*TpYx3CYx3*TrYb3adjYb3/5._dp -       & 
&  3._dp*TpYx3CYx3md2*TrYb3adjYb3/5._dp - 6._dp*TpYx3mHxb32CYx3*TrYb3adjYb3/5._dp -      & 
&  6._dp*TpYx3CYx3*TrYb3ml2adjYb3/5._dp - 6._dp*md2YdadjYd*TrYdadjYd - 2._dp*md2YdadjYd*TrYeadjYe -& 
&  6._dp*md2TpYx3CYx3*TrYuadjYu - 12._dp*TpTYx3CTYx3*TrYuadjYu - 24._dp*mHu2*TpYx3CYx3*TrYuadjYu -& 
&  6._dp*TpYx3CYx3md2*TrYuadjYu - 12._dp*TpYx3mHxb32CYx3*TrYuadjYu - 12._dp*TpYx3CYx3*TrYumq2adjYu -& 
&  3._dp*md2TpYx3CYx3*TrYw3adjYw3 - 6._dp*TpTYx3CTYx3*TrYw3adjYw3 - 12._dp*mHu2*TpYx3CYx3*TrYw3adjYw3 -& 
&  3._dp*TpYx3CYx3md2*TrYw3adjYw3 - 6._dp*TpYx3mHxb32CYx3*TrYw3adjYw3 - 6._dp*TpYx3CYx3*TrYw3ml2adjYw3  
betamd22 =  betamd22- 6._dp*md2TpYx3CYx3*TrYx3adjYx3 - 12._dp*TpTYx3CTYx3*TrYx3adjYx3 - 24._dp*mHu2*TpYx3CYx3*TrYx3adjYx3 -& 
&  6._dp*TpYx3CYx3md2*TrYx3adjYx3 - 12._dp*TpYx3mHxb32CYx3*TrYx3adjYx3 - 12._dp*TpYx3CYx3*TrYx3md2adjYx3 +& 
&  4._dp*g1p2*TYdadjTYd/5._dp + 12._dp*g2p2*TYdadjTYd - 12._dp*TrYdadjYd*TYdadjTYd -     & 
&  4._dp*TrYeadjYe*TYdadjTYd - 12._dp*TrTpYdCTYd*TYdadjYd - 4._dp*TrTpYeCTYe*TYdadjYd -  & 
&  4._dp*TYdadjYdYdadjTYd - 4._dp*TYdadjYuYuadjTYd - 4._dp*TYdTpYdCTYdadjYd -            & 
&  4._dp*TYdTpYuCTYuadjYd - 4._dp*g1p2*MassB*YdadjTYd/5._dp - 12._dp*g2p2*MassWB*YdadjTYd -& 
&  12._dp*TrTYdadjYd*YdadjTYd - 4._dp*TrTYeadjYe*YdadjTYd - 4._dp*YdadjTYdTYdadjYd -     & 
&  4._dp*YdadjTYuTYuadjYd + 4._dp*g1p2*mHd2*YdadjYd/5._dp + 12._dp*g2p2*mHd2*YdadjYd -   & 
&  12._dp*Trmd2YdadjYd*YdadjYd - 4._dp*Trme2YeadjYe*YdadjYd - 12._dp*TrTpTYdCTYd*YdadjYd -& 
&  4._dp*TrTpTYeCTYe*YdadjYd - 24._dp*mHd2*TrYdadjYd*YdadjYd - 12._dp*TrYdmq2adjYd*YdadjYd -& 
&  8._dp*mHd2*TrYeadjYe*YdadjYd - 4._dp*TrYeml2adjYe*YdadjYd + 2._dp*g1p2*YdadjYdmd2/5._dp  
betamd22 =  betamd22+ 6._dp*g2p2*YdadjYdmd2 - 6._dp*TrYdadjYd*YdadjYdmd2 - 2._dp*TrYeadjYe*YdadjYdmd2 -     & 
&  4._dp*YdadjYdmd2YdadjYd - 4._dp*YdadjYdTYdadjTYd - 8._dp*mHd2*YdadjYdYdadjYd -        & 
&  2._dp*YdadjYdYdadjYdmd2 - 4._dp*YdadjYdYdmq2adjYd - 4._dp*YdadjYumu2YuadjYd -         & 
&  4._dp*YdadjYuTYuadjTYd - 4._dp*mHd2*YdadjYuYuadjYd - 4._dp*mHu2*YdadjYuYuadjYd -      & 
&  2._dp*YdadjYuYuadjYdmd2 - 4._dp*YdadjYuYumq2adjYd + 4._dp*g1p2*Ydmq2adjYd/5._dp +     & 
&  12._dp*g2p2*Ydmq2adjYd - 12._dp*TrYdadjYd*Ydmq2adjYd - 4._dp*TrYeadjYe*Ydmq2adjYd -   & 
&  4._dp*Ydmq2adjYdYdadjYd - 4._dp*Ydmq2adjYuYuadjYd + 808._dp*g1p4*id3R*MassB*Conjg(MassB)/75._dp +& 
&  128._dp*g1p2*g3p2*id3R*MassB*Conjg(MassB)/45._dp + 64._dp*g1p2*g3p2*id3R*MassG*Conjg(MassB)/45._dp -& 
&  4._dp*g1p2*TpTYx3CYx3*Conjg(MassB) + 8._dp*g1p2*MassB*TpYx3CYx3*Conjg(MassB) -        & 
&  4._dp*g1p2*TYdadjYd*Conjg(MassB)/5._dp + 8._dp*g1p2*MassB*YdadjYd*Conjg(MassB)/5._dp +& 
&  64._dp*g1p2*g3p2*id3R*MassB*Conjg(MassG)/45._dp + 128._dp*g1p2*g3p2*id3R*MassG*Conjg(MassG)/45._dp  
betamd22 =  betamd22- 128._dp*g3p4*id3R*MassG*Conjg(MassG)/3._dp - 12._dp*g2p2*TpTYx3CYx3*Conjg(MassWB) +   & 
&  24._dp*g2p2*MassWB*TpYx3CYx3*Conjg(MassWB) - 12._dp*g2p2*TYdadjYd*Conjg(MassWB) +     & 
&  24._dp*g2p2*MassWB*YdadjYd*Conjg(MassWB) + 96._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 +& 
&  4._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 +& 
&  4._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3 +& 
&  8._dp*g1p4*id3R*Tr2(1)/15._dp + 32._dp*g3p4*id3R*Tr2(3)/3._dp + 8._dp*g1p2*id3R*Tr3(1)/3._dp

 
Dmd2 = oo16pi2*( betamd21 + oo16pi2 * betamd22 ) 

 
Else 
Dmd2 = oo16pi2* betamd21 
End If 
 
 
Forall(i1=1:3) Dmd2(i1,i1) =  Real(Dmd2(i1,i1),dp) 
!-------------------- 
! mu2 
!-------------------- 
 
betamu21  = 2._dp*mu2YuadjYu + 4._dp*TYuadjTYu + 4._dp*mHu2*YuadjYu + 2._dp*YuadjYumu2 +& 
&  4._dp*Yumq2adjYu - 32._dp*g1p2*id3R*MassB*Conjg(MassB)/15._dp - 32._dp*g3p2*id3R*MassG*Conjg(MassG)& 
& /3._dp - 4._dp*g1p2*id3R*Tr1(1)/3._dp

 
 
If (TwoLoopRGE) Then 
betamu22 = -2._dp*mu2YuadjYdYdadjYu - 2._dp*g1p2*mu2YuadjYu/5._dp + 6._dp*g2p2*mu2YuadjYu -      & 
&  2._dp*mu2YuadjYuYuadjYu - 3._dp*mu2YuadjYu*TrYb3adjYb3/5._dp - 6._dp*mu2YuadjYu*TrYuadjYu -& 
&  3._dp*mu2YuadjYu*TrYw3adjYw3 - 6._dp*mu2YuadjYu*TrYx3adjYx3 - 4._dp*g1p2*TYuadjTYu/5._dp +& 
&  12._dp*g2p2*TYuadjTYu - 6._dp*TrYb3adjYb3*TYuadjTYu/5._dp - 12._dp*TrYuadjYu*TYuadjTYu -& 
&  6._dp*TrYw3adjYw3*TYuadjTYu - 12._dp*TrYx3adjYx3*TYuadjTYu - 4._dp*TYuadjYdYdadjTYu - & 
&  6._dp*TrTpYb3CTYb3*TYuadjYu/5._dp - 12._dp*TrTpYuCTYu*TYuadjYu - 6._dp*TrTpYw3CTYw3*TYuadjYu -& 
&  12._dp*TrTpYx3CTYx3*TYuadjYu - 4._dp*TYuadjYuYuadjTYu - 4._dp*TYuTpYdCTYdadjYu -      & 
&  4._dp*TYuTpYuCTYuadjYu - 4._dp*YuadjTYdTYdadjYu + 4._dp*g1p2*MassB*YuadjTYu/5._dp -   & 
&  12._dp*g2p2*MassWB*YuadjTYu - 6._dp*TrTYb3adjYb3*YuadjTYu/5._dp - 12._dp*TrTYuadjYu*YuadjTYu -& 
&  6._dp*TrTYw3adjYw3*YuadjTYu - 12._dp*TrTYx3adjYx3*YuadjTYu - 4._dp*YuadjTYuTYuadjYu - & 
&  4._dp*YuadjYdmd2YdadjYu - 4._dp*YuadjYdTYdadjTYu - 4._dp*mHd2*YuadjYdYdadjYu  
betamu22 =  betamu22- 4._dp*mHu2*YuadjYdYdadjYu - 2._dp*YuadjYdYdadjYumu2 - 4._dp*YuadjYdYdmq2adjYu -       & 
&  4._dp*g1p2*mHu2*YuadjYu/5._dp + 12._dp*g2p2*mHu2*YuadjYu - 6._dp*TrmHb32Yb3adjYb3*YuadjYu/5._dp -& 
&  6._dp*TrmHw32Yw3adjYw3*YuadjYu - 12._dp*TrmHxb32Yx3adjYx3*YuadjYu - 12._dp*Trmu2YuadjYu*YuadjYu -& 
&  6._dp*TrTpTYb3CTYb3*YuadjYu/5._dp - 12._dp*TrTpTYuCTYu*YuadjYu - 6._dp*TrTpTYw3CTYw3*YuadjYu -& 
&  12._dp*TrTpTYx3CTYx3*YuadjYu - 12._dp*mHu2*TrYb3adjYb3*YuadjYu/5._dp - 6._dp*TrYb3ml2adjYb3*YuadjYu/5._dp -& 
&  24._dp*mHu2*TrYuadjYu*YuadjYu - 12._dp*TrYumq2adjYu*YuadjYu - 12._dp*mHu2*TrYw3adjYw3*YuadjYu -& 
&  6._dp*TrYw3ml2adjYw3*YuadjYu - 24._dp*mHu2*TrYx3adjYx3*YuadjYu - 12._dp*TrYx3md2adjYx3*YuadjYu -& 
&  2._dp*g1p2*YuadjYumu2/5._dp + 6._dp*g2p2*YuadjYumu2 - 3._dp*TrYb3adjYb3*YuadjYumu2/5._dp -& 
&  6._dp*TrYuadjYu*YuadjYumu2 - 3._dp*TrYw3adjYw3*YuadjYumu2 - 6._dp*TrYx3adjYx3*YuadjYumu2 -& 
&  4._dp*YuadjYumu2YuadjYu - 4._dp*YuadjYuTYuadjTYu - 8._dp*mHu2*YuadjYuYuadjYu -        & 
&  2._dp*YuadjYuYuadjYumu2 - 4._dp*YuadjYuYumq2adjYu - 4._dp*Yumq2adjYdYdadjYu  
betamu22 =  betamu22- 4._dp*g1p2*Yumq2adjYu/5._dp + 12._dp*g2p2*Yumq2adjYu - 6._dp*TrYb3adjYb3*Yumq2adjYu/5._dp -& 
&  12._dp*TrYuadjYu*Yumq2adjYu - 6._dp*TrYw3adjYw3*Yumq2adjYu - 12._dp*TrYx3adjYx3*Yumq2adjYu -& 
&  4._dp*Yumq2adjYuYuadjYu + 3424._dp*g1p4*id3R*MassB*Conjg(MassB)/75._dp +              & 
&  512._dp*g1p2*g3p2*id3R*MassB*Conjg(MassB)/45._dp + 256._dp*g1p2*g3p2*id3R*MassG*Conjg(MassB)/45._dp +& 
&  4._dp*g1p2*TYuadjYu*Conjg(MassB)/5._dp - 8._dp*g1p2*MassB*YuadjYu*Conjg(MassB)/5._dp +& 
&  256._dp*g1p2*g3p2*id3R*MassB*Conjg(MassG)/45._dp + 512._dp*g1p2*g3p2*id3R*MassG*Conjg(MassG)/45._dp -& 
&  128._dp*g3p4*id3R*MassG*Conjg(MassG)/3._dp - 12._dp*g2p2*TYuadjYu*Conjg(MassWB) +     & 
&  24._dp*g2p2*MassWB*YuadjYu*Conjg(MassWB) + 96._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 +& 
&  16._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 +& 
&  16._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3 +& 
&  32._dp*g1p4*id3R*Tr2(1)/15._dp + 32._dp*g3p4*id3R*Tr2(3)/3._dp - 16._dp*g1p2*id3R*Tr3(1)/3._dp

 
Dmu2 = oo16pi2*( betamu21 + oo16pi2 * betamu22 ) 

 
Else 
Dmu2 = oo16pi2* betamu21 
End If 
 
 
Forall(i1=1:3) Dmu2(i1,i1) =  Real(Dmu2(i1,i1),dp) 
!-------------------- 
! me2 
!-------------------- 
 
betame21  = 2._dp*me2YeadjYe + 4._dp*TYeadjTYe + 4._dp*mHd2*YeadjYe + 2._dp*YeadjYeme2 +& 
&  4._dp*Yeml2adjYe - 24._dp*g1p2*id3R*MassB*Conjg(MassB)/5._dp + 2._dp*g1p2*id3R*Tr1(1)

 
 
If (TwoLoopRGE) Then 
betame22 = -3._dp*me2YeadjYb3Yb3adjYe/5._dp - 6._dp*g1p2*me2YeadjYe/5._dp + 6._dp*g2p2*me2YeadjYe -& 
&  2._dp*me2YeadjYeYeadjYe - 3._dp*me2YeadjYw3Yw3adjYe - 6._dp*me2YeadjYe*TrYdadjYd -    & 
&  2._dp*me2YeadjYe*TrYeadjYe - 12._dp*g1p2*TYeadjTYe/5._dp + 12._dp*g2p2*TYeadjTYe -    & 
&  12._dp*TrYdadjYd*TYeadjTYe - 4._dp*TrYeadjYe*TYeadjTYe - 6._dp*TYeadjYb3Yb3adjTYe/5._dp -& 
&  12._dp*TrTpYdCTYd*TYeadjYe - 4._dp*TrTpYeCTYe*TYeadjYe - 4._dp*TYeadjYeYeadjTYe -     & 
&  6._dp*TYeadjYw3Yw3adjTYe - 6._dp*TYeTpYb3CTYb3adjYe/5._dp - 4._dp*TYeTpYeCTYeadjYe -  & 
&  6._dp*TYeTpYw3CTYw3adjYe - 6._dp*YeadjTYb3TYb3adjYe/5._dp + 12._dp*g1p2*MassB*YeadjTYe/5._dp -& 
&  12._dp*g2p2*MassWB*YeadjTYe - 12._dp*TrTYdadjYd*YeadjTYe - 4._dp*TrTYeadjYe*YeadjTYe -& 
&  4._dp*YeadjTYeTYeadjYe - 6._dp*YeadjTYw3TYw3adjYe - 6._dp*YeadjYb3mHb32Yb3adjYe/5._dp -& 
&  6._dp*YeadjYb3TYb3adjTYe/5._dp - 6._dp*mHd2*YeadjYb3Yb3adjYe/5._dp - 6._dp*mHu2*YeadjYb3Yb3adjYe/5._dp -& 
&  3._dp*YeadjYb3Yb3adjYeme2/5._dp - 6._dp*YeadjYb3Yb3ml2adjYe/5._dp - 12._dp*g1p2*mHd2*YeadjYe/5._dp  
betame22 =  betame22+ 12._dp*g2p2*mHd2*YeadjYe - 12._dp*Trmd2YdadjYd*YeadjYe - 4._dp*Trme2YeadjYe*YeadjYe - & 
&  12._dp*TrTpTYdCTYd*YeadjYe - 4._dp*TrTpTYeCTYe*YeadjYe - 24._dp*mHd2*TrYdadjYd*YeadjYe -& 
&  12._dp*TrYdmq2adjYd*YeadjYe - 8._dp*mHd2*TrYeadjYe*YeadjYe - 4._dp*TrYeml2adjYe*YeadjYe -& 
&  6._dp*g1p2*YeadjYeme2/5._dp + 6._dp*g2p2*YeadjYeme2 - 6._dp*TrYdadjYd*YeadjYeme2 -    & 
&  2._dp*TrYeadjYe*YeadjYeme2 - 4._dp*YeadjYeme2YeadjYe - 4._dp*YeadjYeTYeadjTYe -       & 
&  8._dp*mHd2*YeadjYeYeadjYe - 2._dp*YeadjYeYeadjYeme2 - 4._dp*YeadjYeYeml2adjYe -       & 
&  6._dp*YeadjYw3mHw32Yw3adjYe - 6._dp*YeadjYw3TYw3adjTYe - 6._dp*mHd2*YeadjYw3Yw3adjYe -& 
&  6._dp*mHu2*YeadjYw3Yw3adjYe - 3._dp*YeadjYw3Yw3adjYeme2 - 6._dp*YeadjYw3Yw3ml2adjYe - & 
&  6._dp*Yeml2adjYb3Yb3adjYe/5._dp - 12._dp*g1p2*Yeml2adjYe/5._dp + 12._dp*g2p2*Yeml2adjYe -& 
&  12._dp*TrYdadjYd*Yeml2adjYe - 4._dp*TrYeadjYe*Yeml2adjYe - 4._dp*Yeml2adjYeYeadjYe -  & 
&  6._dp*Yeml2adjYw3Yw3adjYe + 2808._dp*g1p4*id3R*MassB*Conjg(MassB)/25._dp  
betame22 =  betame22+ 12._dp*g1p2*TYeadjYe*Conjg(MassB)/5._dp - 24._dp*g1p2*MassB*YeadjYe*Conjg(MassB)/5._dp -& 
&  12._dp*g2p2*TYeadjYe*Conjg(MassWB) + 24._dp*g2p2*MassWB*YeadjYe*Conjg(MassWB) +       & 
&  36._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 36._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 +& 
&  24._dp*g1p4*id3R*Tr2(1)/5._dp + 8._dp*g1p2*id3R*Tr3(1)

 
Dme2 = oo16pi2*( betame21 + oo16pi2 * betame22 ) 

 
Else 
Dme2 = oo16pi2* betame21 
End If 
 
 
Forall(i1=1:3) Dme2(i1,i1) =  Real(Dme2(i1,i1),dp) 
!-------------------- 
! mHw32 
!-------------------- 
 
betamHw321  = mHw32Yw3adjYw3 + 2._dp*TYw3adjTYw3 + 2._dp*mHu2*Yw3adjYw3 +             & 
&  Yw3adjYw3mHw32 + 2._dp*Yw3ml2adjYw3 - 16._dp*g2p2*id3R*MassWB*Conjg(MassWB)

 
 
If (TwoLoopRGE) Then 
betamHw322 = -3._dp*mHw32Yw3adjYb3Yb3adjYw3/10._dp - mHw32Yw3adjYeYeadjYw3 + 3._dp*g1p2*mHw32Yw3adjYw3/5._dp -& 
&  g2p2*mHw32Yw3adjYw3 - 3._dp*mHw32Yw3adjYw3Yw3adjYw3/2._dp - 3._dp*mHw32Yw3adjYw3*TrYb3adjYb3/10._dp -& 
&  3._dp*mHw32Yw3adjYw3*TrYuadjYu - 3._dp*mHw32Yw3adjYw3*TrYw3adjYw3/2._dp -             & 
&  3._dp*mHw32Yw3adjYw3*TrYx3adjYx3 + 6._dp*g1p2*TYw3adjTYw3/5._dp - 2._dp*g2p2*TYw3adjTYw3 -& 
&  3._dp*TrYb3adjYb3*TYw3adjTYw3/5._dp - 6._dp*TrYuadjYu*TYw3adjTYw3 - 3._dp*TrYw3adjYw3*TYw3adjTYw3 -& 
&  6._dp*TrYx3adjYx3*TYw3adjTYw3 - 3._dp*TYw3adjYb3Yb3adjTYw3/5._dp - 2._dp*TYw3adjYeYeadjTYw3 -& 
&  3._dp*TrTpYb3CTYb3*TYw3adjYw3/5._dp - 6._dp*TrTpYuCTYu*TYw3adjYw3 - 3._dp*TrTpYw3CTYw3*TYw3adjYw3 -& 
&  6._dp*TrTpYx3CTYx3*TYw3adjYw3 - 3._dp*TYw3adjYw3Yw3adjTYw3 - 3._dp*TYw3TpYb3CTYb3adjYw3/5._dp -& 
&  2._dp*TYw3TpYeCTYeadjYw3 - 3._dp*TYw3TpYw3CTYw3adjYw3 - 3._dp*Yw3adjTYb3TYb3adjYw3/5._dp -& 
&  2._dp*Yw3adjTYeTYeadjYw3 - 6._dp*g1p2*MassB*Yw3adjTYw3/5._dp + 2._dp*g2p2*MassWB*Yw3adjTYw3 -& 
&  3._dp*TrTYb3adjYb3*Yw3adjTYw3/5._dp - 6._dp*TrTYuadjYu*Yw3adjTYw3 - 3._dp*TrTYw3adjYw3*Yw3adjTYw3  
betamHw322 =  betamHw322- 6._dp*TrTYx3adjYx3*Yw3adjTYw3 - 3._dp*Yw3adjTYw3TYw3adjYw3 - 3._dp*Yw3adjYb3mHb32Yb3adjYw3/5._dp -& 
&  3._dp*Yw3adjYb3TYb3adjTYw3/5._dp - 6._dp*mHu2*Yw3adjYb3Yb3adjYw3/5._dp -              & 
&  3._dp*Yw3adjYb3Yb3adjYw3mHw32/10._dp - 3._dp*Yw3adjYb3Yb3ml2adjYw3/5._dp -            & 
&  2._dp*Yw3adjYeme2YeadjYw3 - 2._dp*Yw3adjYeTYeadjTYw3 - 2._dp*mHd2*Yw3adjYeYeadjYw3 -  & 
&  2._dp*mHu2*Yw3adjYeYeadjYw3 - Yw3adjYeYeadjYw3mHw32 - 2._dp*Yw3adjYeYeml2adjYw3 +     & 
&  6._dp*g1p2*mHu2*Yw3adjYw3/5._dp - 2._dp*g2p2*mHu2*Yw3adjYw3 - 3._dp*TrmHb32Yb3adjYb3*Yw3adjYw3/5._dp -& 
&  3._dp*TrmHw32Yw3adjYw3*Yw3adjYw3 - 6._dp*TrmHxb32Yx3adjYx3*Yw3adjYw3 - 6._dp*Trmu2YuadjYu*Yw3adjYw3 -& 
&  3._dp*TrTpTYb3CTYb3*Yw3adjYw3/5._dp - 6._dp*TrTpTYuCTYu*Yw3adjYw3 - 3._dp*TrTpTYw3CTYw3*Yw3adjYw3 -& 
&  6._dp*TrTpTYx3CTYx3*Yw3adjYw3 - 6._dp*mHu2*TrYb3adjYb3*Yw3adjYw3/5._dp -              & 
&  3._dp*TrYb3ml2adjYb3*Yw3adjYw3/5._dp - 12._dp*mHu2*TrYuadjYu*Yw3adjYw3 -              & 
&  6._dp*TrYumq2adjYu*Yw3adjYw3 - 6._dp*mHu2*TrYw3adjYw3*Yw3adjYw3 - 3._dp*TrYw3ml2adjYw3*Yw3adjYw3  
betamHw322 =  betamHw322- 12._dp*mHu2*TrYx3adjYx3*Yw3adjYw3 - 6._dp*TrYx3md2adjYx3*Yw3adjYw3 + 3._dp*g1p2*Yw3adjYw3mHw32/5._dp -& 
&  g2p2*Yw3adjYw3mHw32 - 3._dp*TrYb3adjYb3*Yw3adjYw3mHw32/10._dp - 3._dp*TrYuadjYu*Yw3adjYw3mHw32 -& 
&  3._dp*TrYw3adjYw3*Yw3adjYw3mHw32/2._dp - 3._dp*TrYx3adjYx3*Yw3adjYw3mHw32 -           & 
&  3._dp*Yw3adjYw3mHw32Yw3adjYw3 - 3._dp*Yw3adjYw3TYw3adjTYw3 - 6._dp*mHu2*Yw3adjYw3Yw3adjYw3 -& 
&  3._dp*Yw3adjYw3Yw3adjYw3mHw32/2._dp - 3._dp*Yw3adjYw3Yw3ml2adjYw3 - 3._dp*Yw3ml2adjYb3Yb3adjYw3/5._dp -& 
&  2._dp*Yw3ml2adjYeYeadjYw3 + 6._dp*g1p2*Yw3ml2adjYw3/5._dp - 2._dp*g2p2*Yw3ml2adjYw3 - & 
&  3._dp*TrYb3adjYb3*Yw3ml2adjYw3/5._dp - 6._dp*TrYuadjYu*Yw3ml2adjYw3 - 3._dp*TrYw3adjYw3*Yw3ml2adjYw3 -& 
&  6._dp*TrYx3adjYx3*Yw3ml2adjYw3 - 3._dp*Yw3ml2adjYw3Yw3adjYw3 - 6._dp*g1p2*TYw3adjYw3*Conjg(MassB)/5._dp +& 
&  12._dp*g1p2*MassB*Yw3adjYw3*Conjg(MassB)/5._dp + 208._dp*g2p4*id3R*MassWB*Conjg(MassWB) +& 
&  2._dp*g2p2*TYw3adjYw3*Conjg(MassWB) - 4._dp*g2p2*MassWB*Yw3adjYw3*Conjg(MassWB) +     & 
&  96._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHw3 + 72._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHx3  
betamHw322 =  betamHw322+ 72._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHxb3 + 16._dp*g2p4*id3R*Tr2(2)

 
DmHw32 = oo16pi2*( betamHw321 + oo16pi2 * betamHw322 ) 

 
Else 
DmHw32 = oo16pi2* betamHw321 
End If 
 
 
Forall(i1=1:3) DmHw32(i1,i1) =  Real(DmHw32(i1,i1),dp) 
!-------------------- 
! mHg32 
!-------------------- 
 
betamHg321  = -24._dp*g3p2*id3R*MassG*Conjg(MassG)

 
 
If (TwoLoopRGE) Then 
betamHg322 = 144._dp*g3p4*id3R*MassG*Conjg(MassG) + 216._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 +& 
&  72._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 + 72._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3 +& 
&  24._dp*g3p4*id3R*Tr2(3)

 
DmHg32 = oo16pi2*( betamHg321 + oo16pi2 * betamHg322 ) 

 
Else 
DmHg32 = oo16pi2* betamHg321 
End If 
 
 
Forall(i1=1:3) DmHg32(i1,i1) =  Real(DmHg32(i1,i1),dp) 
!-------------------- 
! mHb32 
!-------------------- 
 
betamHb321  = 3._dp*mHb32Yb3adjYb3/5._dp + 6._dp*TYb3adjTYb3/5._dp + 6._dp*mHu2*Yb3adjYb3/5._dp +& 
&  3._dp*Yb3adjYb3mHb32/5._dp + 6._dp*Yb3ml2adjYb3/5._dp

 
 
If (TwoLoopRGE) Then 
betamHb322 = 9._dp*g1p2*mHb32Yb3adjYb3/25._dp + 9._dp*g2p2*mHb32Yb3adjYb3/5._dp - 9._dp*mHb32Yb3adjYb3Yb3adjYb3/50._dp -& 
&  3._dp*mHb32Yb3adjYeYeadjYb3/5._dp - 9._dp*mHb32Yb3adjYw3Yw3adjYb3/10._dp -            & 
&  9._dp*mHb32Yb3adjYb3*TrYb3adjYb3/50._dp - 9._dp*mHb32Yb3adjYb3*TrYuadjYu/5._dp -      & 
&  9._dp*mHb32Yb3adjYb3*TrYw3adjYw3/10._dp - 9._dp*mHb32Yb3adjYb3*TrYx3adjYx3/5._dp +    & 
&  18._dp*g1p2*TYb3adjTYb3/25._dp + 18._dp*g2p2*TYb3adjTYb3/5._dp - 9._dp*TrYb3adjYb3*TYb3adjTYb3/25._dp -& 
&  18._dp*TrYuadjYu*TYb3adjTYb3/5._dp - 9._dp*TrYw3adjYw3*TYb3adjTYb3/5._dp -            & 
&  18._dp*TrYx3adjYx3*TYb3adjTYb3/5._dp - 9._dp*TrTpYb3CTYb3*TYb3adjYb3/25._dp -         & 
&  18._dp*TrTpYuCTYu*TYb3adjYb3/5._dp - 9._dp*TrTpYw3CTYw3*TYb3adjYb3/5._dp -            & 
&  18._dp*TrTpYx3CTYx3*TYb3adjYb3/5._dp - 9._dp*TYb3adjYb3Yb3adjTYb3/25._dp -            & 
&  6._dp*TYb3adjYeYeadjTYb3/5._dp - 9._dp*TYb3adjYw3Yw3adjTYb3/5._dp - 9._dp*TYb3TpYb3CTYb3adjYb3/25._dp -& 
&  6._dp*TYb3TpYeCTYeadjYb3/5._dp - 9._dp*TYb3TpYw3CTYw3adjYb3/5._dp - 18._dp*g1p2*MassB*Yb3adjTYb3/25._dp  
betamHb322 =  betamHb322- 18._dp*g2p2*MassWB*Yb3adjTYb3/5._dp - 9._dp*TrTYb3adjYb3*Yb3adjTYb3/25._dp -          & 
&  18._dp*TrTYuadjYu*Yb3adjTYb3/5._dp - 9._dp*TrTYw3adjYw3*Yb3adjTYb3/5._dp -            & 
&  18._dp*TrTYx3adjYx3*Yb3adjTYb3/5._dp - 9._dp*Yb3adjTYb3TYb3adjYb3/25._dp -            & 
&  6._dp*Yb3adjTYeTYeadjYb3/5._dp - 9._dp*Yb3adjTYw3TYw3adjYb3/5._dp + 18._dp*g1p2*mHu2*Yb3adjYb3/25._dp +& 
&  18._dp*g2p2*mHu2*Yb3adjYb3/5._dp - 9._dp*TrmHb32Yb3adjYb3*Yb3adjYb3/25._dp -          & 
&  9._dp*TrmHw32Yw3adjYw3*Yb3adjYb3/5._dp - 18._dp*TrmHxb32Yx3adjYx3*Yb3adjYb3/5._dp -   & 
&  18._dp*Trmu2YuadjYu*Yb3adjYb3/5._dp - 9._dp*TrTpTYb3CTYb3*Yb3adjYb3/25._dp -          & 
&  18._dp*TrTpTYuCTYu*Yb3adjYb3/5._dp - 9._dp*TrTpTYw3CTYw3*Yb3adjYb3/5._dp -            & 
&  18._dp*TrTpTYx3CTYx3*Yb3adjYb3/5._dp - 18._dp*mHu2*TrYb3adjYb3*Yb3adjYb3/25._dp -     & 
&  9._dp*TrYb3ml2adjYb3*Yb3adjYb3/25._dp - 36._dp*mHu2*TrYuadjYu*Yb3adjYb3/5._dp -       & 
&  18._dp*TrYumq2adjYu*Yb3adjYb3/5._dp - 18._dp*mHu2*TrYw3adjYw3*Yb3adjYb3/5._dp  
betamHb322 =  betamHb322- 9._dp*TrYw3ml2adjYw3*Yb3adjYb3/5._dp - 36._dp*mHu2*TrYx3adjYx3*Yb3adjYb3/5._dp -      & 
&  18._dp*TrYx3md2adjYx3*Yb3adjYb3/5._dp + 9._dp*g1p2*Yb3adjYb3mHb32/25._dp +            & 
&  9._dp*g2p2*Yb3adjYb3mHb32/5._dp - 9._dp*TrYb3adjYb3*Yb3adjYb3mHb32/50._dp -           & 
&  9._dp*TrYuadjYu*Yb3adjYb3mHb32/5._dp - 9._dp*TrYw3adjYw3*Yb3adjYb3mHb32/10._dp -      & 
&  9._dp*TrYx3adjYx3*Yb3adjYb3mHb32/5._dp - 9._dp*Yb3adjYb3mHb32Yb3adjYb3/25._dp -       & 
&  9._dp*Yb3adjYb3TYb3adjTYb3/25._dp - 18._dp*mHu2*Yb3adjYb3Yb3adjYb3/25._dp -           & 
&  9._dp*Yb3adjYb3Yb3adjYb3mHb32/50._dp - 9._dp*Yb3adjYb3Yb3ml2adjYb3/25._dp -           & 
&  6._dp*Yb3adjYeme2YeadjYb3/5._dp - 6._dp*Yb3adjYeTYeadjTYb3/5._dp - 6._dp*mHd2*Yb3adjYeYeadjYb3/5._dp -& 
&  6._dp*mHu2*Yb3adjYeYeadjYb3/5._dp - 3._dp*Yb3adjYeYeadjYb3mHb32/5._dp -               & 
&  6._dp*Yb3adjYeYeml2adjYb3/5._dp - 9._dp*Yb3adjYw3mHw32Yw3adjYb3/5._dp -               & 
&  9._dp*Yb3adjYw3TYw3adjTYb3/5._dp - 18._dp*mHu2*Yb3adjYw3Yw3adjYb3/5._dp  
betamHb322 =  betamHb322- 9._dp*Yb3adjYw3Yw3adjYb3mHb32/10._dp - 9._dp*Yb3adjYw3Yw3ml2adjYb3/5._dp +            & 
&  18._dp*g1p2*Yb3ml2adjYb3/25._dp + 18._dp*g2p2*Yb3ml2adjYb3/5._dp - 9._dp*TrYb3adjYb3*Yb3ml2adjYb3/25._dp -& 
&  18._dp*TrYuadjYu*Yb3ml2adjYb3/5._dp - 9._dp*TrYw3adjYw3*Yb3ml2adjYb3/5._dp -          & 
&  18._dp*TrYx3adjYx3*Yb3ml2adjYb3/5._dp - 9._dp*Yb3ml2adjYb3Yb3adjYb3/25._dp -          & 
&  6._dp*Yb3ml2adjYeYeadjYb3/5._dp - 9._dp*Yb3ml2adjYw3Yw3adjYb3/5._dp - 18._dp*g1p2*TYb3adjYb3*Conjg(MassB)/25._dp +& 
&  36._dp*g1p2*MassB*Yb3adjYb3*Conjg(MassB)/25._dp - 18._dp*g2p2*TYb3adjYb3*Conjg(MassWB)/5._dp +& 
&  36._dp*g2p2*MassWB*Yb3adjYb3*Conjg(MassWB)/5._dp

 
DmHb32 = oo16pi2*( betamHb321 + oo16pi2 * betamHb322 ) 

 
Else 
DmHb32 = oo16pi2* betamHb321 
End If 
 
 
Forall(i1=1:3) DmHb32(i1,i1) =  Real(DmHb32(i1,i1),dp) 
!-------------------- 
! mHx32 
!-------------------- 
 
betamHx321  = -10._dp*g1p2*id3R*MassB*Conjg(MassB)/3._dp - 32._dp*g3p2*id3R*MassG*Conjg(MassG)& 
& /3._dp - 6._dp*g2p2*id3R*MassWB*Conjg(MassWB) + 5._dp*g1p2*id3R*Tr1(1)/3._dp

 
 
If (TwoLoopRGE) Then 
betamHx322 = 223._dp*g1p4*id3R*MassB*Conjg(MassB)/3._dp + 10._dp*g1p2*g2p2*id3R*MassB*Conjg(MassB) +& 
&  160._dp*g1p2*g3p2*id3R*MassB*Conjg(MassB)/9._dp + 80._dp*g1p2*g3p2*id3R*MassG*Conjg(MassB)/9._dp +& 
&  5._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassB) + 80._dp*g1p2*g3p2*id3R*MassB*Conjg(MassG)/9._dp +& 
&  160._dp*g1p2*g3p2*id3R*MassG*Conjg(MassG)/9._dp + 32._dp*g2p2*g3p2*id3R*MassG*Conjg(MassG) -& 
&  128._dp*g3p4*id3R*MassG*Conjg(MassG)/3._dp + 16._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassG) +& 
&  5._dp*g1p2*g2p2*id3R*MassB*Conjg(MassWB) + 16._dp*g2p2*g3p2*id3R*MassG*Conjg(MassWB) +& 
&  10._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassWB) + 33._dp*g2p4*id3R*MassWB*Conjg(MassWB) +  & 
&  32._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassWB) + 96._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 +& 
&  36._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHw3 + 25._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 +& 
&  32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 + 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHx3 +& 
&  25._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3  
betamHx322 =  betamHx322+ 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHxb3 + 10._dp*g1p4*id3R*Tr2(1)/3._dp +& 
&  6._dp*g2p4*id3R*Tr2(2) + 32._dp*g3p4*id3R*Tr2(3)/3._dp + 20._dp*g1p2*id3R*Tr3(1)/3._dp

 
DmHx32 = oo16pi2*( betamHx321 + oo16pi2 * betamHx322 ) 

 
Else 
DmHx32 = oo16pi2* betamHx321 
End If 
 
 
Forall(i1=1:3) DmHx32(i1,i1) =  Real(DmHx32(i1,i1),dp) 
!-------------------- 
! mHxb32 
!-------------------- 
 
betamHxb321  = mHxb32Yx3adjYx3 + 2._dp*TYx3adjTYx3 + 2._dp*mHu2*Yx3adjYx3 +           & 
&  Yx3adjYx3mHxb32 + 2._dp*Yx3md2adjYx3 - 10._dp*g1p2*id3R*MassB*Conjg(MassB)            & 
& /3._dp - 32._dp*g3p2*id3R*MassG*Conjg(MassG)/3._dp - 6._dp*g2p2*id3R*MassWB*Conjg(MassWB)& 
&  - 5._dp*g1p2*id3R*Tr1(1)/3._dp

 
 
If (TwoLoopRGE) Then 
betamHxb322 = -2._dp*g1p2*mHxb32Yx3adjYx3/5._dp - 2._dp*mHxb32Yx3adjYx3Yx3adjYx3 - 2._dp*mHxb32Yx3CYdTpYdadjYx3 -& 
&  3._dp*mHxb32Yx3adjYx3*TrYb3adjYb3/10._dp - 3._dp*mHxb32Yx3adjYx3*TrYuadjYu -          & 
&  3._dp*mHxb32Yx3adjYx3*TrYw3adjYw3/2._dp - 3._dp*mHxb32Yx3adjYx3*TrYx3adjYx3 -         & 
&  4._dp*g1p2*TYx3adjTYx3/5._dp - 3._dp*TrYb3adjYb3*TYx3adjTYx3/5._dp - 6._dp*TrYuadjYu*TYx3adjTYx3 -& 
&  3._dp*TrYw3adjYw3*TYx3adjTYx3 - 6._dp*TrYx3adjYx3*TYx3adjTYx3 - 3._dp*TrTpYb3CTYb3*TYx3adjYx3/5._dp -& 
&  6._dp*TrTpYuCTYu*TYx3adjYx3 - 3._dp*TrTpYw3CTYw3*TYx3adjYx3 - 6._dp*TrTpYx3CTYx3*TYx3adjYx3 -& 
&  4._dp*TYx3adjYx3Yx3adjTYx3 - 4._dp*TYx3CYdTpYdadjTYx3 - 4._dp*TYx3TpYx3CTYx3adjYx3 -  & 
&  4._dp*TYx3YdadjTYdadjYx3 + 4._dp*g1p2*MassB*Yx3adjTYx3/5._dp - 3._dp*TrTYb3adjYb3*Yx3adjTYx3/5._dp -& 
&  6._dp*TrTYuadjYu*Yx3adjTYx3 - 3._dp*TrTYw3adjYw3*Yx3adjTYx3 - 6._dp*TrTYx3adjYx3*Yx3adjTYx3 -& 
&  4._dp*Yx3adjTYx3TYx3adjYx3 - 4._dp*g1p2*mHu2*Yx3adjYx3/5._dp - 3._dp*TrmHb32Yb3adjYb3*Yx3adjYx3/5._dp -& 
&  3._dp*TrmHw32Yw3adjYw3*Yx3adjYx3 - 6._dp*TrmHxb32Yx3adjYx3*Yx3adjYx3 - 6._dp*Trmu2YuadjYu*Yx3adjYx3  
betamHxb322 =  betamHxb322- 3._dp*TrTpTYb3CTYb3*Yx3adjYx3/5._dp - 6._dp*TrTpTYuCTYu*Yx3adjYx3 - 3._dp*TrTpTYw3CTYw3*Yx3adjYx3 -& 
&  6._dp*TrTpTYx3CTYx3*Yx3adjYx3 - 6._dp*mHu2*TrYb3adjYb3*Yx3adjYx3/5._dp -              & 
&  3._dp*TrYb3ml2adjYb3*Yx3adjYx3/5._dp - 12._dp*mHu2*TrYuadjYu*Yx3adjYx3 -              & 
&  6._dp*TrYumq2adjYu*Yx3adjYx3 - 6._dp*mHu2*TrYw3adjYw3*Yx3adjYx3 - 3._dp*TrYw3ml2adjYw3*Yx3adjYx3 -& 
&  12._dp*mHu2*TrYx3adjYx3*Yx3adjYx3 - 6._dp*TrYx3md2adjYx3*Yx3adjYx3 - 2._dp*g1p2*Yx3adjYx3mHxb32/5._dp -& 
&  3._dp*TrYb3adjYb3*Yx3adjYx3mHxb32/10._dp - 3._dp*TrYuadjYu*Yx3adjYx3mHxb32 -          & 
&  3._dp*TrYw3adjYw3*Yx3adjYx3mHxb32/2._dp - 3._dp*TrYx3adjYx3*Yx3adjYx3mHxb32 -         & 
&  4._dp*Yx3adjYx3mHxb32Yx3adjYx3 - 4._dp*Yx3adjYx3TYx3adjTYx3 - 8._dp*mHu2*Yx3adjYx3Yx3adjYx3 -& 
&  2._dp*Yx3adjYx3Yx3adjYx3mHxb32 - 4._dp*Yx3adjYx3Yx3md2adjYx3 - 4._dp*Yx3CTYdTpTYdadjYx3 -& 
&  4._dp*Yx3CYdmq2TpYdadjYx3 - 4._dp*mHd2*Yx3CYdTpYdadjYx3 - 4._dp*mHu2*Yx3CYdTpYdadjYx3 -& 
&  2._dp*Yx3CYdTpYdadjYx3mHxb32 - 4._dp*Yx3CYdTpYdmd2adjYx3 - 4._dp*g1p2*Yx3md2adjYx3/5._dp  
betamHxb322 =  betamHxb322- 3._dp*TrYb3adjYb3*Yx3md2adjYx3/5._dp - 6._dp*TrYuadjYu*Yx3md2adjYx3 - 3._dp*TrYw3adjYw3*Yx3md2adjYx3 -& 
&  6._dp*TrYx3adjYx3*Yx3md2adjYx3 - 4._dp*Yx3md2adjYx3Yx3adjYx3 - 4._dp*Yx3md2CYdTpYdadjYx3 -& 
&  4._dp*Yx3TYdadjYdadjTYx3 + 223._dp*g1p4*id3R*MassB*Conjg(MassB)/3._dp +               & 
&  10._dp*g1p2*g2p2*id3R*MassB*Conjg(MassB) + 160._dp*g1p2*g3p2*id3R*MassB*Conjg(MassB)/9._dp +& 
&  80._dp*g1p2*g3p2*id3R*MassG*Conjg(MassB)/9._dp + 5._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassB) +& 
&  4._dp*g1p2*TYx3adjYx3*Conjg(MassB)/5._dp - 8._dp*g1p2*MassB*Yx3adjYx3*Conjg(MassB)/5._dp +& 
&  80._dp*g1p2*g3p2*id3R*MassB*Conjg(MassG)/9._dp + 160._dp*g1p2*g3p2*id3R*MassG*Conjg(MassG)/9._dp +& 
&  32._dp*g2p2*g3p2*id3R*MassG*Conjg(MassG) - 128._dp*g3p4*id3R*MassG*Conjg(MassG)/3._dp +& 
&  16._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassG) + 5._dp*g1p2*g2p2*id3R*MassB*Conjg(MassWB) +& 
&  16._dp*g2p2*g3p2*id3R*MassG*Conjg(MassWB) + 10._dp*g1p2*g2p2*id3R*MassWB*Conjg(MassWB) +& 
&  33._dp*g2p4*id3R*MassWB*Conjg(MassWB) + 32._dp*g2p2*g3p2*id3R*MassWB*Conjg(MassWB)  
betamHxb322 =  betamHxb322+ 96._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHg3 + 36._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHw3 +& 
&  25._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHx3 + 32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHx3 +& 
&  27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHx3 + 25._dp*g1p4*id3R*MassB*Conjg(MassB)*NGHxb3 +& 
&  32._dp*g3p4*id3R*MassG*Conjg(MassG)*NGHxb3 + 27._dp*g2p4*id3R*MassWB*Conjg(MassWB)*NGHxb3 +& 
&  10._dp*g1p4*id3R*Tr2(1)/3._dp + 6._dp*g2p4*id3R*Tr2(2) + 32._dp*g3p4*id3R*Tr2(3)/3._dp -& 
&  20._dp*g1p2*id3R*Tr3(1)/3._dp

 
DmHxb32 = oo16pi2*( betamHxb321 + oo16pi2 * betamHxb322 ) 

 
Else 
DmHxb32 = oo16pi2* betamHxb321 
End If 
 
 
Forall(i1=1:3) DmHxb32(i1,i1) =  Real(DmHxb32(i1,i1),dp) 
!-------------------- 
! MassB 
!-------------------- 
 
betaMassB1  = 66._dp*g1p2*MassB/5._dp + 5._dp*g1p2*MassB*NGHx3 +       & 
&  5._dp*g1p2*MassB*NGHxb3

 
 
If (TwoLoopRGE) Then 
betaMassB2 = 796._dp*g1p4*MassB/25._dp + 54._dp*g1p2*g2p2*MassB/5._dp + 176._dp*g1p2*g3p2*MassB/5._dp +& 
&  176._dp*g1p2*g3p2*MassG/5._dp + 54._dp*g1p2*g2p2*MassWB/5._dp + 18._dp*g1p2*TrTYb3adjYb3/25._dp +& 
&  28._dp*g1p2*TrTYdadjYd/5._dp + 36._dp*g1p2*TrTYeadjYe/5._dp + 52._dp*g1p2*TrTYuadjYu/5._dp +& 
&  24._dp*g1p2*TrTYw3adjYw3/5._dp + 76._dp*g1p2*TrTYx3adjYx3/5._dp - 18._dp*g1p2*MassB*TrYb3adjYb3/25._dp -& 
&  28._dp*g1p2*MassB*TrYdadjYd/5._dp - 36._dp*g1p2*MassB*TrYeadjYe/5._dp -               & 
&  52._dp*g1p2*MassB*TrYuadjYu/5._dp - 24._dp*g1p2*MassB*TrYw3adjYw3/5._dp -             & 
&  76._dp*g1p2*MassB*TrYx3adjYx3/5._dp + 50._dp*g1p4*MassB*NGHx3/3._dp +  & 
&  15._dp*g1p2*g2p2*MassB*NGHx3 + 80._dp*g1p2*g3p2*MassB*NGHx3/3._dp +& 
&  80._dp*g1p2*g3p2*MassG*NGHx3/3._dp + 15._dp*g1p2*g2p2*MassWB*NGHx3 +& 
&  50._dp*g1p4*MassB*NGHxb3/3._dp + 15._dp*g1p2*g2p2*MassB*NGHxb3 +& 
&  80._dp*g1p2*g3p2*MassB*NGHxb3/3._dp + 80._dp*g1p2*g3p2*MassG*NGHxb3/3._dp  
betaMassB2 =  betaMassB2+ 15._dp*g1p2*g2p2*MassWB*NGHxb3

 
DMassB = oo16pi2*( betaMassB1 + oo16pi2 * betaMassB2 ) 

 
Else 
DMassB = oo16pi2* betaMassB1 
End If 
 
 
!-------------------- 
! MassWB 
!-------------------- 
 
betaMassWB1  = 2._dp*g2p2*MassWB + 4._dp*g2p2*MassWB*NGHw3 +           & 
&  3._dp*g2p2*MassWB*NGHx3 + 3._dp*g2p2*MassWB*NGHxb3

 
 
If (TwoLoopRGE) Then 
betaMassWB2 = 18._dp*g1p2*g2p2*MassB/5._dp + 48._dp*g2p2*g3p2*MassG + 18._dp*g1p2*g2p2*MassWB/5._dp +& 
&  100._dp*g2p4*MassWB + 48._dp*g2p2*g3p2*MassWB + 6._dp*g2p2*TrTYb3adjYb3/5._dp +       & 
&  12._dp*g2p2*TrTYdadjYd + 4._dp*g2p2*TrTYeadjYe + 12._dp*g2p2*TrTYuadjYu +             & 
&  56._dp*g2p2*TrTYw3adjYw3/3._dp + 12._dp*g2p2*TrTYx3adjYx3 - 6._dp*g2p2*MassWB*TrYb3adjYb3/5._dp -& 
&  12._dp*g2p2*MassWB*TrYdadjYd - 4._dp*g2p2*MassWB*TrYeadjYe - 12._dp*g2p2*MassWB*TrYuadjYu -& 
&  56._dp*g2p2*MassWB*TrYw3adjYw3/3._dp - 12._dp*g2p2*MassWB*TrYx3adjYx3 +               & 
&  96._dp*g2p4*MassWB*NGHw3 + 5._dp*g1p2*g2p2*MassB*NGHx3 +& 
&  16._dp*g2p2*g3p2*MassG*NGHx3 + 5._dp*g1p2*g2p2*MassWB*NGHx3 +& 
&  42._dp*g2p4*MassWB*NGHx3 + 16._dp*g2p2*g3p2*MassWB*NGHx3 +& 
&  5._dp*g1p2*g2p2*MassB*NGHxb3 + 16._dp*g2p2*g3p2*MassG*NGHxb3 +& 
&  5._dp*g1p2*g2p2*MassWB*NGHxb3 + 42._dp*g2p4*MassWB*NGHxb3  
betaMassWB2 =  betaMassWB2+ 16._dp*g2p2*g3p2*MassWB*NGHxb3

 
DMassWB = oo16pi2*( betaMassWB1 + oo16pi2 * betaMassWB2 ) 

 
Else 
DMassWB = oo16pi2* betaMassWB1 
End If 
 
 
!-------------------- 
! MassG 
!-------------------- 
 
betaMassG1  = -6._dp*g3p2*MassG + 6._dp*g3p2*MassG*NGHg3 +             & 
&  2._dp*g3p2*MassG*NGHx3 + 2._dp*g3p2*MassG*NGHxb3

 
 
If (TwoLoopRGE) Then 
betaMassG2 = 22._dp*g1p2*g3p2*MassB/5._dp + 22._dp*g1p2*g3p2*MassG/5._dp + 18._dp*g2p2*g3p2*MassG +& 
&  56._dp*g3p4*MassG + 18._dp*g2p2*g3p2*MassWB + 8._dp*g3p2*TrTYdadjYd + 8._dp*g3p2*TrTYuadjYu +& 
&  8._dp*g3p2*TrTYx3adjYx3 - 8._dp*g3p2*MassG*TrYdadjYd - 8._dp*g3p2*MassG*TrYuadjYu -   & 
&  8._dp*g3p2*MassG*TrYx3adjYx3 + 216._dp*g3p4*MassG*NGHg3 +              & 
&  10._dp*g1p2*g3p2*MassB*NGHx3/3._dp + 10._dp*g1p2*g3p2*MassG*NGHx3/3._dp +& 
&  6._dp*g2p2*g3p2*MassG*NGHx3 + 136._dp*g3p4*MassG*NGHx3/3._dp +& 
&  6._dp*g2p2*g3p2*MassWB*NGHx3 + 10._dp*g1p2*g3p2*MassB*NGHxb3/3._dp +& 
&  10._dp*g1p2*g3p2*MassG*NGHxb3/3._dp + 6._dp*g2p2*g3p2*MassG*NGHxb3 +& 
&  136._dp*g3p4*MassG*NGHxb3/3._dp + 6._dp*g2p2*g3p2*MassWB*NGHxb3

 
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
& 3._dp*TrYb3adjYw3Yw3adjYb3/4._dp - 6._dp*TrYdadjYdTpYx3CYx3 - &
& 3._dp*TrYuadjYdYdadjYu + 4._dp*g1p2*TrYuadjYu/5._dp + 16._dp*g3p2*TrYuadjYu &
& - 9._dp*TrYuadjYuYuadjYu - 27._dp*TrYw3adjYb3Yb3adjYw3/40._dp - &
& 3._dp*TrYw3adjYeYeadjYw3/2._dp + 6._dp*g2p2*TrYw3adjYw3 - &
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

Subroutine GToParameters117(g,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM)

Implicit None 
Real(dp), Intent(in) :: g(117) 
Real(dp),Intent(out) :: g1,g2,g3

Complex(dp),Intent(out) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)            & 
& ,L1,L2,MTM

Integer i1, i2, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'GToParameters117' 
 
g1= g(1) 
g2= g(2) 
g3= g(3) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yt(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ys(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
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
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+4) = Real(Yu(i1,i2), dp) 
g(SumI+5) = Aimag(Yu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+22) = Real(Yd(i1,i2), dp) 
g(SumI+23) = Aimag(Yd(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+40) = Real(Ye(i1,i2), dp) 
g(SumI+41) = Aimag(Ye(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+58) = Real(Yt(i1,i2), dp) 
g(SumI+59) = Aimag(Yt(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+76) = Real(Ys(i1,i2), dp) 
g(SumI+77) = Aimag(Ys(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
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

Subroutine rge117(len, T, GY, F) 
Implicit None 
Integer, Intent(in) :: len 
Real(dp), Intent(in) :: T, GY(len) 
Real(dp), Intent(out) :: F(len) 
Integer :: i1,i2 
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

Complex(dp) :: YdadjYu(3,3),YeadjYz(3,3),YeCYt(3,3),YtadjYe(3,3),YtadjYz(3,3)         & 
& ,YuadjYd(3,3),YzadjYe(3,3),YzCYt(3,3),CYeTYz(3,3),CYtTYz(3,3),CYuTYd(3,3)              & 
& ,YeadjYzYd(3,3),YeadjYzYs(3,3),YsCYzYt(3,3),YtadjYzYd(3,3),YtadjYzYs(3,3)              & 
& ,YtCYtTYz(3,3),YuadjYdYs(3,3),YuadjYdYz(3,3),adjYdYdadjYd(3,3),adjYdYdadjYu(3,3)       & 
& ,adjYdYsCYs(3,3),adjYdYzadjYz(3,3),adjYeYeadjYe(3,3),adjYeYeadjYz(3,3),adjYeYeCYt(3,3) & 
& ,adjYuYuadjYd(3,3),adjYuYuadjYu(3,3),adjYzYdadjYd(3,3),adjYzYsCYs(3,3),adjYzYzadjYe(3,3)& 
& ,adjYzYzadjYz(3,3),adjYzYzCYt(3,3),CYsYdadjYd(3,3),CYsYsCYs(3,3),CYsYzadjYz(3,3)       & 
& ,CYtYtadjYe(3,3),CYtYtadjYz(3,3),CYtYtCYt(3,3),TYdCYdTYd(3,3),TYdCYsYd(3,3)            & 
& ,TYdCYsYs(3,3),TYdCYsYz(3,3),TYdCYzYt(3,3),TYeCYeTYz(3,3),TYuCYuTYd(3,3)               & 
& ,TYzCYsYd(3,3),TYzCYsYs(3,3),TYzCYsYz(3,3),TYzCYzTYz(3,3),YdadjYdYdadjYd(3,3)          & 
& ,YdadjYdYsCYs(3,3),YdadjYdYzadjYz(3,3),YdadjYuYuadjYd(3,3),YeadjYeYeadjYe(3,3)         & 
& ,YeadjYzYzadjYe(3,3),YeCYtYtadjYe(3,3),YsCYsYdadjYd(3,3),YsCYsYsCYs(3,3)               & 
& ,YsCYsYzadjYz(3,3),YtadjYeYeCYt(3,3),YtadjYzYzCYt(3,3),YtCYtYtCYt(3,3),YuadjYdYdadjYu(3,3)& 
& ,YuadjYuYuadjYu(3,3),YzadjYeYeadjYz(3,3),YzadjYzYdadjYd(3,3),YzadjYzYsCYs(3,3)         & 
& ,YzadjYzYzadjYz(3,3),YzCYtYtadjYz(3,3),adjYdYdadjYdYd(3,3),adjYdYdadjYdYs(3,3)         & 
& ,adjYdYdadjYdYz(3,3),adjYdYdadjYuYu(3,3),adjYdYsCYsYd(3,3),adjYdYzadjYzYd(3,3)         & 
& ,adjYeYeadjYeYe(3,3),adjYeYeadjYzYd(3,3),adjYeYeadjYzYs(3,3),adjYeYeadjYzYz(3,3)       & 
& ,adjYeYeCYtYt(3,3),adjYuYuadjYdYd(3,3),adjYuYuadjYdYs(3,3),adjYuYuadjYdYz(3,3)         & 
& ,adjYuYuadjYuYu(3,3),adjYzYdadjYdYz(3,3),adjYzYsCYsYz(3,3),adjYzYzadjYeYe(3,3)         & 
& ,adjYzYzadjYzYd(3,3),adjYzYzadjYzYs(3,3),adjYzYzadjYzYz(3,3),adjYzYzCYtYt(3,3)         & 
& ,CYdTYdCYdTYd(3,3),CYdTYdCYsYd(3,3),CYdTYdCYsYs(3,3),CYdTYdCYsYz(3,3),CYdTYdCYzYt(3,3) & 
& ,CYdTYuCYuTYd(3,3),CYeTYeCYeYt(3,3),CYsYdadjYdYs(3,3),CYsYsCYsYd(3,3),CYsYsCYsYs(3,3)  & 
& ,CYsYsCYsYz(3,3),CYsYsCYzYt(3,3),CYsYzadjYzYs(3,3),CYtYtadjYeYe(3,3),CYtYtadjYzYd(3,3) & 
& ,CYtYtadjYzYs(3,3),CYtYtadjYzYz(3,3),CYtYtCYtYt(3,3),CYtTYeCYeYt(3,3),CYtTYzCYzYt(3,3) & 
& ,CYzYtCYtTYz(3,3),CYzTYeCYeTYz(3,3),CYzTYzCYsYd(3,3),CYzTYzCYsYs(3,3),CYzTYzCYsYz(3,3) & 
& ,CYzTYzCYzYt(3,3),CYzTYzCYzTYz(3,3),YdadjYdYdadjYdYd(3,3),YdadjYdYdadjYdYs(3,3)        & 
& ,YdadjYdYdadjYdYz(3,3),YdadjYdYsCYsYd(3,3),YdadjYdYzadjYzYd(3,3),YdadjYuYuadjYdYd(3,3) & 
& ,YdadjYuYuadjYdYs(3,3),YdadjYuYuadjYdYz(3,3),YdadjYuYuadjYuYu(3,3),YeadjYeYeadjYeYe(3,3)& 
& ,YeadjYzYdadjYdYz(3,3),YeadjYzYsCYsYz(3,3),YeadjYzYzadjYeYe(3,3),YeadjYzYzadjYzYz(3,3) & 
& ,YeCYtYtadjYeYe(3,3),YeCYtYtCYtYt(3,3),YeCYtTYeCYeYt(3,3),YeCYtTYzCYzYt(3,3)           & 
& ,YsCYdTYdCYdTYd(3,3),YsCYdTYdCYsYd(3,3),YsCYdTYdCYsYs(3,3),YsCYdTYdCYsYz(3,3)          & 
& ,YsCYdTYuCYuTYd(3,3),YsCYsYdadjYdYs(3,3),YsCYsYsCYsYd(3,3),YsCYsYsCYsYs(3,3)           & 
& ,YsCYsYsCYsYz(3,3),YsCYsYzadjYzYs(3,3),YsCYzYtCYtTYz(3,3),YsCYzTYeCYeTYz(3,3)          & 
& ,YsCYzTYzCYsYd(3,3),YsCYzTYzCYsYs(3,3),YsCYzTYzCYsYz(3,3),YsCYzTYzCYzTYz(3,3)          & 
& ,YtadjYeYeadjYeYe(3,3),YtadjYeYeCYtYt(3,3),YtadjYzYdadjYdYz(3,3),YtadjYzYsCYsYz(3,3)   & 
& ,YtadjYzYzadjYzYz(3,3),YtadjYzYzCYtYt(3,3),YtCYtYtCYtYt(3,3),YtCYtTYeCYeYt(3,3)        & 
& ,YtCYtTYzCYzYt(3,3),YuadjYdYdadjYdYd(3,3),YuadjYdYdadjYuYu(3,3),YuadjYdYsCYsYd(3,3)    & 
& ,YuadjYdYzadjYzYd(3,3),YuadjYuYuadjYuYu(3,3),YzadjYeYeadjYeYe(3,3),YzadjYeYeadjYzYd(3,3)& 
& ,YzadjYeYeadjYzYs(3,3),YzadjYeYeadjYzYz(3,3),YzadjYzYdadjYdYz(3,3),YzadjYzYsCYsYz(3,3) & 
& ,YzadjYzYzadjYzYd(3,3),YzadjYzYzadjYzYs(3,3),YzadjYzYzadjYzYz(3,3),YzCYtYtadjYzYd(3,3) & 
& ,YzCYtYtadjYzYs(3,3),YzCYtYtadjYzYz(3,3),YzCYtYtCYtYt(3,3),YzCYtTYeCYeYt(3,3)          & 
& ,YzCYtTYzCYzYt(3,3),TYeCYeTYeCYeYt(3,3),TYzCYdTYdCYzYt(3,3),TYzCYsYsCYzYt(3,3)         & 
& ,TYzCYzTYzCYzYt(3,3)

Complex(dp) :: TrYdadjYd,TrYeadjYe,TrYsCYs,TrYtCYt,TrYuadjYu,TrYzadjYz,               & 
& TrCYsYs,TrCYtYt

Complex(dp) :: TrYdadjYdYdadjYd,TrYdadjYdYsCYs,TrYdadjYdYzadjYz,TrYdadjYuYuadjYd,     & 
& TrYeadjYeYeadjYe,TrYeadjYzYzadjYe,TrYeCYtYtadjYe,TrYsCYsYdadjYd,TrYsCYsYsCYs,          & 
& TrYsCYsYzadjYz,TrYtadjYeYeCYt,TrYtadjYzYzCYt,TrYtCYtYtCYt,TrYuadjYdYdadjYu,            & 
& TrYuadjYuYuadjYu,TrYzadjYeYeadjYz,TrYzadjYzYdadjYd,TrYzadjYzYsCYs,TrYzadjYzYzadjYz,    & 
& TrYzCYtYtadjYz,TrCYsYdadjYdYs,TrCYsYsCYsYs,TrCYsYzadjYzYs,TrCYtYtadjYeYe,              & 
& TrCYtYtadjYzYz,TrCYtYtCYtYt

 Real(dp) :: g12, g13, g22, g23, g32, g33

Iname = Iname +1 
NameOfUnit(Iname) = 'rge117' 
 
OnlyDiagonal = .Not.GenerationMixing 


q = t 
 
Call GToParameters117(gy,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM)

 g12 = g1**2
 g13 = g1*g12
 g22 = g2**2
 g23 = g2*g22
 g32 = g3**2
 g33 = g3*g32

Call Adjungate(Yu,adjYu)
Call Adjungate(Yd,adjYd)
Call Adjungate(Ye,adjYe)
Call Adjungate(Yt,adjYt)
Call Adjungate(Ys,adjYs)
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
 TrCYsYs = Real(cTrace(CYsYs),dp) 
 TrCYtYt = Real(cTrace(CYtYt),dp) 


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
 YtadjYeYeCYt = Matmul2(Yt,adjYeYeCYt,OnlyDiagonal) 
 YtadjYzYzCYt = Matmul2(Yt,adjYzYzCYt,OnlyDiagonal) 
 YtCYtYtCYt = Matmul2(Yt,CYtYtCYt,OnlyDiagonal) 
 YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
 YzadjYeYeadjYz = Matmul2(Yz,adjYeYeadjYz,OnlyDiagonal) 
 YzadjYzYdadjYd = Matmul2(Yz,adjYzYdadjYd,OnlyDiagonal) 
 YzadjYzYsCYs = Matmul2(Yz,adjYzYsCYs,OnlyDiagonal) 
 YzadjYzYzadjYz = Matmul2(Yz,adjYzYzadjYz,OnlyDiagonal) 
 YzCYtYtadjYz = Matmul2(Yz,CYtYtadjYz,OnlyDiagonal) 
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
 TrYdadjYdYzadjYz = cTrace(YdadjYdYzadjYz) 
 TrYdadjYuYuadjYd = cTrace(YdadjYuYuadjYd) 
 TrYeadjYeYeadjYe = cTrace(YeadjYeYeadjYe) 
 TrYeadjYzYzadjYe = cTrace(YeadjYzYzadjYe) 
 TrYeCYtYtadjYe = cTrace(YeCYtYtadjYe) 
 TrYsCYsYdadjYd = cTrace(YsCYsYdadjYd) 
 TrYsCYsYsCYs = cTrace(YsCYsYsCYs) 
 TrYsCYsYzadjYz = cTrace(YsCYsYzadjYz) 
 TrYtadjYeYeCYt = cTrace(YtadjYeYeCYt) 
 TrYtadjYzYzCYt = cTrace(YtadjYzYzCYt) 
 TrYtCYtYtCYt = cTrace(YtCYtYtCYt) 
 TrYuadjYdYdadjYu = cTrace(YuadjYdYdadjYu) 
 TrYuadjYuYuadjYu = cTrace(YuadjYuYuadjYu) 
 TrYzadjYeYeadjYz = cTrace(YzadjYeYeadjYz) 
 TrYzadjYzYdadjYd = cTrace(YzadjYzYdadjYd) 
 TrYzadjYzYsCYs = cTrace(YzadjYzYsCYs) 
 TrYzadjYzYzadjYz = cTrace(YzadjYzYzadjYz) 
 TrYzCYtYtadjYz = cTrace(YzCYtYtadjYz) 
 TrCYsYdadjYdYs = cTrace(CYsYdadjYdYs) 
 TrCYsYsCYsYs = cTrace(CYsYsCYsYs) 
 TrCYsYzadjYzYs = cTrace(CYsYzadjYzYs) 
 TrCYtYtadjYeYe = cTrace(CYtYtadjYeYe) 
 TrCYtYtadjYzYz = cTrace(CYtYtadjYzYz) 
 TrCYtYtCYtYt = cTrace(CYtYtCYtYt) 
End If 
 
 
!-------------------- 
! g1 
!-------------------- 
 
betag11 = 68._dp*g13/5._dp

 
If (TwoLoopRGE) Then 
betag12 = 1502._dp*g1**5/75._dp + 174._dp*g13*g22/5._dp + 184._dp*g13*g32/3._dp -& 
&  21._dp*g13*TrCYsYs/5._dp - 9._dp*g13*TrCYtYt/2._dp - 14._dp*g13*TrYdadjYd/5._dp -& 
&  18._dp*g13*TrYeadjYe/5._dp - 3._dp*g13*TrYsCYs/5._dp - 9._dp*g13*TrYtCYt/10._dp -& 
&  26._dp*g13*TrYuadjYu/5._dp - 14._dp*g13*TrYzadjYz/5._dp - 27._dp*g13*L1*Conjg(L1)& 
& /5._dp - 27._dp*g13*L2*Conjg(L2)/5._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21 = 8._dp*g23

 
If (TwoLoopRGE) Then 
betag22 = 58._dp*g12*g23/5._dp + 94._dp*g2**5 + 40._dp*g23*g32 -              & 
&  35._dp*g23*TrCYtYt/6._dp - 6._dp*g23*TrYdadjYd - 2._dp*g23*TrYeadjYe -          & 
&  7._dp*g23*TrYtCYt/6._dp - 6._dp*g23*TrYuadjYu - 6._dp*g23*TrYzadjYz -           & 
&  7._dp*g23*L1*Conjg(L1) - 7._dp*g23*L2*Conjg(L2)

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31 = 4._dp*g33

 
If (TwoLoopRGE) Then 
betag32 = 23._dp*g12*g33/3._dp + 15._dp*g22*g33 + 400._dp*g3**5/3._dp -       & 
&  63._dp*g33*TrCYsYs/8._dp - 4._dp*g33*TrYdadjYd - 9._dp*g33*TrYsCYs/8._dp -      & 
&  4._dp*g33*TrYuadjYu - 4._dp*g33*TrYzadjYz

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYu1(i1,i2) = -13._dp*g12*Yu(i1,i2)/15._dp - 3._dp*g22*Yu(i1,i2)               & 
&  - 16._dp*g32*Yu(i1,i2)/3._dp + 3._dp*TrYuadjYu*Yu(i1,i2) + 3._dp*L2*Conjg(L2)       & 
& *Yu(i1,i2) + YuadjYdYd(i1,i2) + 3._dp*YuadjYuYu(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYu2(i1,i2) = 5473._dp*g1**4*Yu(i1,i2)/450._dp + g12*g22*Yu(i1,i2)             & 
&  + 57._dp*g2**4*Yu(i1,i2)/2._dp + 136._dp*g12*g32*Yu(i1,i2)/45._dp +               & 
&  8._dp*g22*g32*Yu(i1,i2) + 320._dp*g3**4*Yu(i1,i2)/9._dp - 3._dp*TrYuadjYdYdadjYu*Yu(i1,i2)& 
&  + 4._dp*g12*TrYuadjYu*Yu(i1,i2)/5._dp + 16._dp*g32*TrYuadjYu*Yu(i1,i2)            & 
&  - 9._dp*TrYuadjYuYuadjYu*Yu(i1,i2) + 18._dp*g12*L2*Conjg(L2)*Yu(i1,i2)              & 
& /5._dp + 12._dp*g22*L2*Conjg(L2)*Yu(i1,i2) - 9._dp*L2*TrYuadjYu*Conjg(L2)            & 
& *Yu(i1,i2) - 12._dp*L2**2*Conjg(L2)**2*Yu(i1,i2) + 2._dp*g12*YuadjYdYd(i1,i2)        & 
& /5._dp - 3._dp*TrYdadjYd*YuadjYdYd(i1,i2) - TrYeadjYe*YuadjYdYd(i1,i2) -               & 
&  3._dp*L1*Conjg(L1)*YuadjYdYd(i1,i2) - 2._dp*YuadjYdYdadjYdYd(i1,i2) - 2._dp*YuadjYdYdadjYuYu(i1,i2)& 
&  - 4._dp*YuadjYdYsCYsYd(i1,i2) - 2._dp*YuadjYdYzadjYzYd(i1,i2) + 2._dp*g12*YuadjYuYu(i1,i2)& 
& /5._dp + 6._dp*g22*YuadjYuYu(i1,i2) - 9._dp*TrYuadjYu*YuadjYuYu(i1,i2)               & 
&  - 9._dp*L2*Conjg(L2)*YuadjYuYu(i1,i2) - 4._dp*YuadjYuYuadjYuYu(i1,i2)
 End Do 
End Do 

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYd1(i1,i2) = -7._dp*g12*Yd(i1,i2)/15._dp - 3._dp*g22*Yd(i1,i2) -              & 
&  16._dp*g32*Yd(i1,i2)/3._dp + 3._dp*TrYdadjYd*Yd(i1,i2) + TrYeadjYe*Yd(i1,i2)        & 
&  + 3._dp*L1*Conjg(L1)*Yd(i1,i2) + 3._dp*YdadjYdYd(i1,i2) + YdadjYuYu(i1,i2)            & 
&  + 4._dp*YsCYsYd(i1,i2) + 2._dp*YzadjYzYd(i1,i2)
 End Do 
End Do 



 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYd2(i1,i2) = 581._dp*g1**4*Yd(i1,i2)/90._dp + g12*g22*Yd(i1,i2)               & 
&  + 57._dp*g2**4*Yd(i1,i2)/2._dp + 8._dp*g12*g32*Yd(i1,i2)/9._dp + 8._dp*g22*g32*Yd(i1,i2)& 
&  + 320._dp*g3**4*Yd(i1,i2)/9._dp - 2._dp*g12*TrYdadjYd*Yd(i1,i2)/5._dp +             & 
&  16._dp*g32*TrYdadjYd*Yd(i1,i2) - 9._dp*TrYdadjYdYdadjYd*Yd(i1,i2) - 3._dp*TrYdadjYuYuadjYd*Yd(i1,i2)& 
&  + 6._dp*g12*TrYeadjYe*Yd(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*Yd(i1,i2)            & 
&  - 3._dp*TrYeadjYzYzadjYe*Yd(i1,i2) - 3._dp*TrYeCYtYtadjYe*Yd(i1,i2) - 12._dp*TrYsCYsYdadjYd*Yd(i1,i2)& 
&  - 6._dp*TrYzadjYzYdadjYd*Yd(i1,i2) + 18._dp*g12*L1*Conjg(L1)*Yd(i1,i2)              & 
& /5._dp + 12._dp*g22*L1*Conjg(L1)*Yd(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)              & 
& *Yd(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)*Yd(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)& 
& *Yd(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)*Yd(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)       & 
& **2*Yd(i1,i2) + 4._dp*g12*YdadjYdYd(i1,i2)/5._dp + 6._dp*g22*YdadjYdYd(i1,i2)      & 
&  - 9._dp*TrYdadjYd*YdadjYdYd(i1,i2) - 3._dp*TrYeadjYe*YdadjYdYd(i1,i2) -               & 
&  9._dp*L1*Conjg(L1)*YdadjYdYd(i1,i2) - 4._dp*YdadjYdYdadjYdYd(i1,i2) - 4._dp*YdadjYdYsCYsYd(i1,i2)& 
&  - 2._dp*YdadjYdYzadjYzYd(i1,i2) + 4._dp*g12*YdadjYuYu(i1,i2)/5._dp - 3._dp*TrYuadjYu*YdadjYuYu(i1,i2)& 
&  - 3._dp*L2*Conjg(L2)*YdadjYuYu(i1,i2) - 2._dp*YdadjYuYuadjYdYd(i1,i2) -               & 
&  2._dp*YdadjYuYuadjYuYu(i1,i2) - 8._dp*YsCYdTYdCYsYd(i1,i2) + 32._dp*g12*YsCYsYd(i1,i2)& 
& /15._dp + 80._dp*g32*YsCYsYd(i1,i2)/3._dp - 3._dp*TrCYsYs*YsCYsYd(i1,i2)             & 
&  - TrYsCYs*YsCYsYd(i1,i2) - 16._dp*YsCYsYsCYsYd(i1,i2) - 8._dp*YsCYzTYzCYsYd(i1,i2)    & 
&  - 2._dp*YzadjYeYeadjYzYd(i1,i2) + 2._dp*g12*YzadjYzYd(i1,i2)/5._dp + 6._dp*g22*YzadjYzYd(i1,i2)& 
&  - 2._dp*TrYzadjYz*YzadjYzYd(i1,i2) - 6._dp*YzadjYzYzadjYzYd(i1,i2) - 6._dp*YzCYtYtadjYzYd(i1,i2)
 End Do 
End Do 

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYe1(i1,i2) = -9._dp*g12*Ye(i1,i2)/5._dp - 3._dp*g22*Ye(i1,i2) +               & 
&  3._dp*TrYdadjYd*Ye(i1,i2) + TrYeadjYe*Ye(i1,i2) + 3._dp*L1*Conjg(L1)*Ye(i1,i2)        & 
&  + 3._dp*YeadjYeYe(i1,i2) + 3._dp*YeadjYzYz(i1,i2) + 3._dp*YeCYtYt(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYe2(i1,i2) = 261._dp*g1**4*Ye(i1,i2)/10._dp + 9._dp*g12*g22*Ye(i1,i2)         & 
& /5._dp + 57._dp*g2**4*Ye(i1,i2)/2._dp - 2._dp*g12*TrYdadjYd*Ye(i1,i2)/5._dp +        & 
&  16._dp*g32*TrYdadjYd*Ye(i1,i2) - 9._dp*TrYdadjYdYdadjYd*Ye(i1,i2) - 3._dp*TrYdadjYuYuadjYd*Ye(i1,i2)& 
&  + 6._dp*g12*TrYeadjYe*Ye(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*Ye(i1,i2)            & 
&  - 3._dp*TrYeadjYzYzadjYe*Ye(i1,i2) - 3._dp*TrYeCYtYtadjYe*Ye(i1,i2) - 12._dp*TrYsCYsYdadjYd*Ye(i1,i2)& 
&  - 6._dp*TrYzadjYzYdadjYd*Ye(i1,i2) + 18._dp*g12*L1*Conjg(L1)*Ye(i1,i2)              & 
& /5._dp + 12._dp*g22*L1*Conjg(L1)*Ye(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)              & 
& *Ye(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)*Ye(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)& 
& *Ye(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)*Ye(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)       & 
& **2*Ye(i1,i2) + 6._dp*g22*YeadjYeYe(i1,i2) - 9._dp*TrYdadjYd*YeadjYeYe(i1,i2)        & 
&  - 3._dp*TrYeadjYe*YeadjYeYe(i1,i2) - 9._dp*L1*Conjg(L1)*YeadjYeYe(i1,i2)              & 
&  - 4._dp*YeadjYeYeadjYeYe(i1,i2) - 6._dp*YeadjYzYdadjYdYz(i1,i2) - 12._dp*YeadjYzYsCYsYz(i1,i2)& 
&  - 2._dp*g12*YeadjYzYz(i1,i2)/5._dp + 16._dp*g32*YeadjYzYz(i1,i2) - 3._dp*TrYzadjYz*YeadjYzYz(i1,i2)& 
&  - 6._dp*YeadjYzYzadjYeYe(i1,i2) - 6._dp*YeadjYzYzadjYzYz(i1,i2) - 3._dp*YeCYtTYeCYeYt(i1,i2)& 
&  - 9._dp*YeCYtTYzCYzYt(i1,i2) + 18._dp*g12*YeCYtYt(i1,i2)/5._dp + 12._dp*g22*YeCYtYt(i1,i2)& 
&  - 9._dp*TrCYtYt*YeCYtYt(i1,i2)/4._dp - 3._dp*TrYtCYt*YeCYtYt(i1,i2)/4._dp -           & 
&  3._dp*L1*Conjg(L1)*YeCYtYt(i1,i2) - 6._dp*YeCYtYtadjYeYe(i1,i2) - 9._dp*YeCYtYtCYtYt(i1,i2)
 End Do 
End Do 

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yt 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYt1(i1,i2) = (sqrt2*TYeCYeYt(i1,i2) + 3._dp*sqrt2*TYzCYzYt(i1,i2)     & 
&  - 9._dp*sqrt2*g12*Yt(i1,i2)/5._dp - 7._dp*sqrt2*g22*Yt(i1,i2)         & 
&  + 3._dp*TrCYtYt*Yt(i1,i2)/(2._dp*sqrt2) + (TrYtCYt*Yt(i1,i2))/(2._dp*sqrt2)& 
&  + sqrt2*L1*Conjg(L1)*Yt(i1,i2) + sqrt2*YtadjYeYe(i1,i2) + 3._dp*sqrt2& 
& *YtadjYzYz(i1,i2) + 6._dp*sqrt2*YtCYtYt(i1,i2))/sqrt2
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYt2(i1,i2) = (-2._dp*sqrt2*TYeCYeTYeCYeYt(i1,i2) + 6._dp*sqrt2        & 
& *g12*TYeCYeYt(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*TYeCYeYt(i1,i2)             & 
&  - sqrt2*TrYeadjYe*TYeCYeYt(i1,i2) - 3._dp*sqrt2*L1*Conjg(L1)              & 
& *TYeCYeYt(i1,i2) - 6._dp*sqrt2*TYzCYdTYdCYzYt(i1,i2) - 12._dp*sqrt2        & 
& *TYzCYsYsCYzYt(i1,i2) - 6._dp*sqrt2*TYzCYzTYzCYzYt(i1,i2) - 2._dp*sqrt2    & 
& *g12*TYzCYzYt(i1,i2)/5._dp + 16._dp*sqrt2*g32*TYzCYzYt(i1,i2) -              & 
&  3._dp*sqrt2*TrYzadjYz*TYzCYzYt(i1,i2) + 417._dp*g1**4*Yt(i1,i2)/(25._dp*sqrt2)& 
&  + 444._dp*sqrt2*g1**4*Yt(i1,i2)/25._dp + 57._dp*sqrt2*g12*g22*Yt(i1,i2)& 
& /5._dp + 57._dp*(g2**4*Yt(i1,i2))/sqrt2 + 48._dp*sqrt2*g2**4*Yt(i1,i2)     & 
&  - 9._dp*g12*TrCYtYt*Yt(i1,i2)/(10._dp*sqrt2) - 3._dp*g22*TrCYtYt*Yt(i1,i2)  & 
& /(2._dp*sqrt2) - sqrt2*TrCYtYtadjYeYe*Yt(i1,i2) - 3._dp*sqrt2        & 
& *TrCYtYtadjYzYz*Yt(i1,i2) - 19._dp*TrCYtYtCYtYt*Yt(i1,i2)/(2._dp*sqrt2)          & 
&  - (TrYeCYtYtadjYe*Yt(i1,i2))/sqrt2 - (TrYtadjYeYeCYt*Yt(i1,i2))/sqrt2     & 
&  - 3._dp*(TrYtadjYzYzCYt*Yt(i1,i2))/sqrt2 - 3._dp*g12*TrYtCYt*Yt(i1,i2)        & 
& /(10._dp*sqrt2) - (g22*TrYtCYt*Yt(i1,i2))/(2._dp*sqrt2) - 5._dp*TrYtCYtYtCYt*Yt(i1,i2)& 
& /(2._dp*sqrt2) - 3._dp*(TrYzCYtYtadjYz*Yt(i1,i2))/sqrt2 - 3._dp*sqrt2& 
& *g12*L1*Conjg(L1)*Yt(i1,i2)/5._dp - sqrt2*g22*L1*Conjg(L1)*Yt(i1,i2)         & 
&  - 6._dp*sqrt2*L1*TrYdadjYd*Conjg(L1)*Yt(i1,i2) - 2._dp*sqrt2              & 
& *L1*TrYeadjYe*Conjg(L1)*Yt(i1,i2) - 6._dp*sqrt2*L1**2*Conjg(L1)**2*Yt(i1,i2)     & 
&  + 6._dp*sqrt2*g12*YtadjYeYe(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*YtadjYeYe(i1,i2)& 
&  - sqrt2*TrYeadjYe*YtadjYeYe(i1,i2) - 3._dp*sqrt2*L1*Conjg(L1)             & 
& *YtadjYeYe(i1,i2) - 2._dp*sqrt2*YtadjYeYeadjYeYe(i1,i2) - 3._dp*sqrt2      & 
& *YtadjYeYeCYtYt(i1,i2) - 6._dp*sqrt2*YtadjYzYdadjYdYz(i1,i2) - 12._dp*sqrt2& 
& *YtadjYzYsCYsYz(i1,i2) - 2._dp*sqrt2*g12*YtadjYzYz(i1,i2)/5._dp +              & 
&  16._dp*sqrt2*g32*YtadjYzYz(i1,i2) - 3._dp*sqrt2*TrYzadjYz*YtadjYzYz(i1,i2)& 
&  - 6._dp*sqrt2*YtadjYzYzadjYzYz(i1,i2) - 9._dp*sqrt2*YtadjYzYzCYtYt(i1,i2) & 
&  - 3._dp*sqrt2*YtCYtTYeCYeYt(i1,i2) - 9._dp*sqrt2*YtCYtTYzCYzYt(i1,i2)     & 
&  + 36._dp*sqrt2*g12*YtCYtYt(i1,i2)/5._dp + 24._dp*sqrt2*g22*YtCYtYt(i1,i2)& 
&  - 9._dp*(TrCYtYt*YtCYtYt(i1,i2))/sqrt2 - 3._dp*(TrYtCYt*YtCYtYt(i1,i2))         & 
& /sqrt2 - 6._dp*sqrt2*L1*Conjg(L1)*YtCYtYt(i1,i2) - 18._dp*sqrt2      & 
& *YtCYtYtCYtYt(i1,i2))/sqrt2
 End Do 
End Do 

 
DYt = oo16pi2*( betaYt1 + oo16pi2 * betaYt2 ) 

 
Else 
DYt = oo16pi2* betaYt1 
End If 
 
 
!-------------------- 
! Ys 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYs1(i1,i2) = (2._dp*sqrt2*YdadjYdYs(i1,i2) - 4._dp*sqrt2              & 
& *g12*Ys(i1,i2)/5._dp - 12._dp*sqrt2*g32*Ys(i1,i2) + 3._dp*TrCYsYs*Ys(i1,i2)  & 
& /(2._dp*sqrt2) + (TrYsCYs*Ys(i1,i2))/(2._dp*sqrt2) + 2._dp*sqrt2     & 
& *YsCYdTYd(i1,i2) + 8._dp*sqrt2*YsCYsYs(i1,i2) + 2._dp*sqrt2*YsCYzTYz(i1,i2)& 
&  + 2._dp*sqrt2*YzadjYzYs(i1,i2))/sqrt2
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYs2(i1,i2) = (-2._dp*sqrt2*YdadjYdYdadjYdYs(i1,i2) + 2._dp*sqrt2      & 
& *g12*YdadjYdYs(i1,i2)/5._dp + 6._dp*sqrt2*g22*YdadjYdYs(i1,i2)               & 
&  - 6._dp*sqrt2*TrYdadjYd*YdadjYdYs(i1,i2) - 2._dp*sqrt2*TrYeadjYe*YdadjYdYs(i1,i2)& 
&  - 6._dp*sqrt2*L1*Conjg(L1)*YdadjYdYs(i1,i2) - 2._dp*sqrt2*YdadjYuYuadjYdYs(i1,i2)& 
&  + 56._dp*sqrt2*g1**4*Ys(i1,i2)/5._dp + 128._dp*sqrt2*g12*g32*Ys(i1,i2)& 
& /15._dp + 320._dp*sqrt2*g3**4*Ys(i1,i2)/3._dp - 2._dp*sqrt2*TrCYsYdadjYdYs*Ys(i1,i2)& 
&  - (sqrt2*g12*TrCYsYs*Ys(i1,i2))/5._dp - sqrt2*g32*TrCYsYs*Ys(i1,i2)   & 
&  - 7._dp*(TrCYsYsCYsYs*Ys(i1,i2))/sqrt2 - 3._dp*sqrt2*TrCYsYsCYsYs*Ys(i1,i2)& 
&  - 2._dp*sqrt2*TrCYsYzadjYzYs*Ys(i1,i2) - sqrt2*TrYdadjYdYsCYs*Ys(i1,i2)   & 
&  - (sqrt2*g12*TrYsCYs*Ys(i1,i2))/15._dp - (sqrt2*g32*TrYsCYs*Ys(i1,i2))& 
& /3._dp - sqrt2*TrYsCYsYdadjYd*Ys(i1,i2) - (TrYsCYsYsCYs*Ys(i1,i2))               & 
& /sqrt2 - sqrt2*TrYsCYsYsCYs*Ys(i1,i2) - sqrt2*TrYsCYsYzadjYz*Ys(i1,i2)& 
&  - sqrt2*TrYzadjYzYsCYs*Ys(i1,i2) + 2._dp*sqrt2*g12*YsCYdTYd(i1,i2)      & 
& /5._dp + 6._dp*sqrt2*g22*YsCYdTYd(i1,i2) - 6._dp*sqrt2*TrYdadjYd*YsCYdTYd(i1,i2)& 
&  - 2._dp*sqrt2*TrYeadjYe*YsCYdTYd(i1,i2) - 6._dp*sqrt2*L1*Conjg(L1)        & 
& *YsCYdTYd(i1,i2) - 2._dp*sqrt2*YsCYdTYdCYdTYd(i1,i2) - 8._dp*sqrt2         & 
& *YsCYdTYdCYsYs(i1,i2) - 2._dp*sqrt2*YsCYdTYuCYuTYd(i1,i2) - 8._dp*sqrt2    & 
& *YsCYsYdadjYdYs(i1,i2) + 64._dp*sqrt2*g12*YsCYsYs(i1,i2)/15._dp +              & 
&  160._dp*sqrt2*g32*YsCYsYs(i1,i2)/3._dp - 6._dp*sqrt2*TrCYsYs*YsCYsYs(i1,i2)& 
&  - 2._dp*sqrt2*TrYsCYs*YsCYsYs(i1,i2) - 32._dp*sqrt2*YsCYsYsCYsYs(i1,i2)   & 
&  - 8._dp*sqrt2*YsCYsYzadjYzYs(i1,i2) - 2._dp*sqrt2*YsCYzTYeCYeTYz(i1,i2)   & 
&  + 2._dp*sqrt2*g12*YsCYzTYz(i1,i2)/5._dp + 6._dp*sqrt2*g22*YsCYzTYz(i1,i2)& 
&  - 2._dp*sqrt2*TrYzadjYz*YsCYzTYz(i1,i2) - 8._dp*sqrt2*YsCYzTYzCYsYs(i1,i2)& 
&  - 6._dp*sqrt2*YsCYzTYzCYzTYz(i1,i2) - 6._dp*sqrt2*YsCYzYtCYtTYz(i1,i2)    & 
&  - 2._dp*sqrt2*YzadjYeYeadjYzYs(i1,i2) + 2._dp*sqrt2*g12*YzadjYzYs(i1,i2)& 
& /5._dp + 6._dp*sqrt2*g22*YzadjYzYs(i1,i2) - 2._dp*sqrt2*TrYzadjYz*YzadjYzYs(i1,i2)& 
&  - 6._dp*sqrt2*YzadjYzYzadjYzYs(i1,i2) - 6._dp*sqrt2*YzCYtYtadjYzYs(i1,i2))& 
& /sqrt2
 End Do 
End Do 

 
DYs = oo16pi2*( betaYs1 + oo16pi2 * betaYs2 ) 

 
Else 
DYs = oo16pi2* betaYs1 
End If 
 
 
!-------------------- 
! Yz 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYz1(i1,i2) = 2._dp*YdadjYdYz(i1,i2) + 4._dp*YsCYsYz(i1,i2) - 7._dp*g12*Yz(i1,i2)& 
& /15._dp - 3._dp*g22*Yz(i1,i2) - 16._dp*g32*Yz(i1,i2)/3._dp + TrYzadjYz*Yz(i1,i2)   & 
&  + YzadjYeYe(i1,i2) + 5._dp*YzadjYzYz(i1,i2) + 3._dp*YzCYtYt(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYz2(i1,i2) = -2._dp*YdadjYdYdadjYdYz(i1,i2) + 2._dp*g12*YdadjYdYz(i1,i2)        & 
& /5._dp + 6._dp*g22*YdadjYdYz(i1,i2) - 6._dp*TrYdadjYd*YdadjYdYz(i1,i2)               & 
&  - 2._dp*TrYeadjYe*YdadjYdYz(i1,i2) - 6._dp*L1*Conjg(L1)*YdadjYdYz(i1,i2)              & 
&  - 2._dp*YdadjYuYuadjYdYz(i1,i2) - 8._dp*YsCYdTYdCYsYz(i1,i2) - 16._dp*YsCYsYsCYsYz(i1,i2)& 
&  + 32._dp*g12*YsCYsYz(i1,i2)/15._dp + 80._dp*g32*YsCYsYz(i1,i2)/3._dp -            & 
&  3._dp*TrCYsYs*YsCYsYz(i1,i2) - TrYsCYs*YsCYsYz(i1,i2) - 8._dp*YsCYzTYzCYsYz(i1,i2)    & 
&  + 581._dp*g1**4*Yz(i1,i2)/90._dp + g12*g22*Yz(i1,i2) + 57._dp*g2**4*Yz(i1,i2)     & 
& /2._dp + 8._dp*g12*g32*Yz(i1,i2)/9._dp + 8._dp*g22*g32*Yz(i1,i2)               & 
&  + 320._dp*g3**4*Yz(i1,i2)/9._dp - 2._dp*TrYdadjYdYzadjYz*Yz(i1,i2) - 4._dp*TrYsCYsYzadjYz*Yz(i1,i2)& 
&  - TrYzadjYeYeadjYz*Yz(i1,i2) + 2._dp*g12*TrYzadjYz*Yz(i1,i2)/5._dp - 5._dp*TrYzadjYzYzadjYz*Yz(i1,i2)& 
&  - 3._dp*TrYzCYtYtadjYz*Yz(i1,i2) + 6._dp*g12*YzadjYeYe(i1,i2)/5._dp -               & 
&  3._dp*TrYdadjYd*YzadjYeYe(i1,i2) - TrYeadjYe*YzadjYeYe(i1,i2) - 3._dp*L1*Conjg(L1)    & 
& *YzadjYeYe(i1,i2) - 2._dp*YzadjYeYeadjYeYe(i1,i2) - 2._dp*YzadjYeYeadjYzYz(i1,i2)      & 
&  - 6._dp*YzadjYzYdadjYdYz(i1,i2) - 12._dp*YzadjYzYsCYsYz(i1,i2) + 6._dp*g22*YzadjYzYz(i1,i2)& 
&  + 16._dp*g32*YzadjYzYz(i1,i2) - 5._dp*TrYzadjYz*YzadjYzYz(i1,i2) - 12._dp*YzadjYzYzadjYzYz(i1,i2)& 
&  - 3._dp*YzCYtTYeCYeYt(i1,i2) - 9._dp*YzCYtTYzCYzYt(i1,i2) + 18._dp*g12*YzCYtYt(i1,i2)& 
& /5._dp + 12._dp*g22*YzCYtYt(i1,i2) - 9._dp*TrCYtYt*YzCYtYt(i1,i2)/4._dp -            & 
&  3._dp*TrYtCYt*YzCYtYt(i1,i2)/4._dp - 3._dp*L1*Conjg(L1)*YzCYtYt(i1,i2) -              & 
&  6._dp*YzCYtYtadjYzYz(i1,i2) - 9._dp*YzCYtYtCYtYt(i1,i2)
 End Do 
End Do 

 
DYz = oo16pi2*( betaYz1 + oo16pi2 * betaYz2 ) 

 
Else 
DYz = oo16pi2* betaYz1 
End If 
 
 
!-------------------- 
! L1 
!-------------------- 
 
betaL11 = (-9._dp*sqrt2*g12*L1/5._dp - 7._dp*sqrt2*g22*L1 + 3._dp*L1*TrCYtYt/(2._dp*sqrt2)& 
&  + 6._dp*sqrt2*L1*TrYdadjYd + 2._dp*sqrt2*L1*TrYeadjYe + (L1*TrYtCYt)      & 
& /(2._dp*sqrt2) + 7._dp*sqrt2*L1**2*Conjg(L1))/sqrt2

 
If (TwoLoopRGE) Then 
betaL12 = (417._dp*g1**4*L1/(25._dp*sqrt2) + 444._dp*sqrt2*g1**4*L1/25._dp + 57._dp*sqrt2& 
& *g12*g22*L1/5._dp + 57._dp*(g2**4*L1)/sqrt2 + 48._dp*sqrt2             & 
& *g2**4*L1 - 9._dp*g12*L1*TrCYtYt/(10._dp*sqrt2) - 3._dp*g22*L1*TrCYtYt/(2._dp*sqrt2)& 
&  - sqrt2*L1*TrCYtYtadjYeYe - 3._dp*sqrt2*L1*TrCYtYtadjYzYz -               & 
&  19._dp*L1*TrCYtYtCYtYt/(2._dp*sqrt2) - 4._dp*sqrt2*g12*L1*TrYdadjYd/5._dp +& 
&  32._dp*sqrt2*g32*L1*TrYdadjYd - 18._dp*sqrt2*L1*TrYdadjYdYdadjYd -      & 
&  6._dp*sqrt2*L1*TrYdadjYuYuadjYd + 12._dp*sqrt2*g12*L1*TrYeadjYe/5._dp - & 
&  6._dp*sqrt2*L1*TrYeadjYeYeadjYe - 6._dp*sqrt2*L1*TrYeadjYzYzadjYe -       & 
&  7._dp*(L1*TrYeCYtYtadjYe)/sqrt2 - 3._dp*sqrt2*L1*TrYeCYtYtadjYe -         & 
&  24._dp*sqrt2*L1*TrYsCYsYdadjYd - (L1*TrYtadjYeYeCYt)/sqrt2 -              & 
&  3._dp*(L1*TrYtadjYzYzCYt)/sqrt2 - 3._dp*g12*L1*TrYtCYt/(10._dp*sqrt2)   & 
&  - (g22*L1*TrYtCYt)/(2._dp*sqrt2) - 5._dp*L1*TrYtCYtYtCYt/(2._dp*sqrt2)  & 
&  - 12._dp*sqrt2*L1*TrYzadjYzYdadjYd - 3._dp*(L1*TrYzCYtYtadjYz)/sqrt2      & 
&  + 33._dp*sqrt2*g12*L1**2*Conjg(L1)/5._dp + 23._dp*sqrt2*g22*L1**2*Conjg(L1)& 
&  - 9._dp*(L1**2*TrCYtYt*Conjg(L1))/sqrt2 - 24._dp*sqrt2*L1**2*TrYdadjYd*Conjg(L1)& 
&  - 8._dp*sqrt2*L1**2*TrYeadjYe*Conjg(L1) - 3._dp*(L1**2*TrYtCYt*Conjg(L1))       & 
& /sqrt2 - 30._dp*sqrt2*L1**3*Conjg(L1)**2)/sqrt2

 
DL1 = oo16pi2*( betaL11 + oo16pi2 * betaL12 ) 

 
Else 
DL1 = oo16pi2* betaL11 
End If 
 
 
!-------------------- 
! L2 
!-------------------- 
 
betaL21 = (-9._dp*sqrt2*g12*L2/5._dp - 7._dp*sqrt2*g22*L2 + 6._dp*sqrt2& 
& *L2*TrYuadjYu + 7._dp*sqrt2*L2**2*Conjg(L2))/sqrt2

 
If (TwoLoopRGE) Then 
betaL22 = (417._dp*g1**4*L2/(25._dp*sqrt2) + 444._dp*sqrt2*g1**4*L2/25._dp + 57._dp*sqrt2& 
& *g12*g22*L2/5._dp + 57._dp*(g2**4*L2)/sqrt2 + 48._dp*sqrt2             & 
& *g2**4*L2 - 6._dp*sqrt2*L2*TrYuadjYdYdadjYu + 8._dp*sqrt2*g12*L2*TrYuadjYu/5._dp +& 
&  32._dp*sqrt2*g32*L2*TrYuadjYu - 18._dp*sqrt2*L2*TrYuadjYuYuadjYu +      & 
&  33._dp*sqrt2*g12*L2**2*Conjg(L2)/5._dp + 23._dp*sqrt2*g22*L2**2*Conjg(L2)& 
&  - 24._dp*sqrt2*L2**2*TrYuadjYu*Conjg(L2) - 30._dp*sqrt2*L2**3*Conjg(L2)   & 
& **2)/sqrt2

 
DL2 = oo16pi2*( betaL21 + oo16pi2 * betaL22 ) 

 
Else 
DL2 = oo16pi2* betaL21 
End If 
 
 
!-------------------- 
! MTM 
!-------------------- 
 
betaMTM1 = -12._dp*g12*MTM/5._dp - 8._dp*g22*MTM + 3._dp*MTM*TrCYtYt/4._dp +      & 
&  (MTM*TrYtCYt)/4._dp + L1*MTM*Conjg(L1) + L2*MTM*Conjg(L2)

 
If (TwoLoopRGE) Then 
betaMTM2 = 888._dp*g1**4*MTM/25._dp + 96._dp*g12*g22*MTM/5._dp + 96._dp*g2**4*MTM -& 
&  9._dp*g12*MTM*TrCYtYt/20._dp - 3._dp*g22*MTM*TrCYtYt/4._dp - MTM*TrCYtYtadjYeYe - & 
&  3._dp*MTM*TrCYtYtadjYzYz - 19._dp*MTM*TrCYtYtCYtYt/4._dp - (MTM*TrYeCYtYtadjYe)       & 
& /2._dp - (MTM*TrYtadjYeYeCYt)/2._dp - 3._dp*MTM*TrYtadjYzYzCYt/2._dp - 3._dp*g12*MTM*TrYtCYt/20._dp -& 
&  (g22*MTM*TrYtCYt)/4._dp - 5._dp*MTM*TrYtCYtYtCYt/4._dp - 3._dp*MTM*TrYzCYtYtadjYz/2._dp -& 
&  3._dp*g12*L1*MTM*Conjg(L1)/5._dp - g22*L1*MTM*Conjg(L1) - 6._dp*L1*MTM*TrYdadjYd*Conjg(L1)& 
&  - 2._dp*L1*MTM*TrYeadjYe*Conjg(L1) - 6._dp*L1**2*MTM*Conjg(L1)**2 - 3._dp*g12*L2*MTM*Conjg(L2)& 
& /5._dp - g22*L2*MTM*Conjg(L2) - 6._dp*L2*MTM*TrYuadjYu*Conjg(L2) - 6._dp*L2**2*MTM*Conjg(L2)**2

 
DMTM = oo16pi2*( betaMTM1 + oo16pi2 * betaMTM2 ) 

 
Else 
DMTM = oo16pi2* betaMTM1 
End If 
 

 
Call ParametersToG117(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYt,DYs,DYz,DL1,DL2,DMTM,f)

Iname = Iname - 1 
 
End Subroutine rge117

Subroutine GToParameters353(g,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,               & 
& MZM,MSM,AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2,              & 
& mHu2,md2,mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG)

Implicit None 
Real(dp), Intent(in) :: g(353) 
Real(dp),Intent(out) :: g1,g2,g3,mHd2,mHu2

Complex(dp),Intent(out) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)            & 
& ,L1,L2,MTM,mue,MZM,MSM,AYu(3,3),AYd(3,3),AYe(3,3),AYt(3,3),AYs(3,3),AYz(3,3)           & 
& ,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3)              & 
& ,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG

Integer i1, i2, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'GToParameters353' 
 
g1= g(1) 
g2= g(2) 
g3= g(3) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yu(i1,i2) = Cmplx( g(SumI+4), g(SumI+5), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yd(i1,i2) = Cmplx( g(SumI+22), g(SumI+23), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ye(i1,i2) = Cmplx( g(SumI+40), g(SumI+41), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yt(i1,i2) = Cmplx( g(SumI+58), g(SumI+59), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Ys(i1,i2) = Cmplx( g(SumI+76), g(SumI+77), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
Yz(i1,i2) = Cmplx( g(SumI+94), g(SumI+95), dp) 
End Do 
 End Do 
 
L1= Cmplx(g(112),g(113),dp) 
L2= Cmplx(g(114),g(115),dp) 
MTM= Cmplx(g(116),g(117),dp) 
mue= Cmplx(g(118),g(119),dp) 
MZM= Cmplx(g(120),g(121),dp) 
MSM= Cmplx(g(122),g(123),dp) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYu(i1,i2) = Cmplx( g(SumI+124), g(SumI+125), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYd(i1,i2) = Cmplx( g(SumI+142), g(SumI+143), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYe(i1,i2) = Cmplx( g(SumI+160), g(SumI+161), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYt(i1,i2) = Cmplx( g(SumI+178), g(SumI+179), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYs(i1,i2) = Cmplx( g(SumI+196), g(SumI+197), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
AYz(i1,i2) = Cmplx( g(SumI+214), g(SumI+215), dp) 
End Do 
 End Do 
 
AL1= Cmplx(g(232),g(233),dp) 
AL2= Cmplx(g(234),g(235),dp) 
Amue= Cmplx(g(236),g(237),dp) 
AMTM= Cmplx(g(238),g(239),dp) 
AMZM= Cmplx(g(240),g(241),dp) 
AMSM= Cmplx(g(242),g(243),dp) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mq2(i1,i2) = Cmplx( g(SumI+244), g(SumI+245), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
ml2(i1,i2) = Cmplx( g(SumI+262), g(SumI+263), dp) 
End Do 
 End Do 
 
mHd2= g(280) 
mHu2= g(281) 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
md2(i1,i2) = Cmplx( g(SumI+282), g(SumI+283), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
mu2(i1,i2) = Cmplx( g(SumI+300), g(SumI+301), dp) 
End Do 
 End Do 
 
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
me2(i1,i2) = Cmplx( g(SumI+318), g(SumI+319), dp) 
End Do 
 End Do 
 
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

Subroutine ParametersToG353(g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,MZM,             & 
& MSM,AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2,mHu2,             & 
& md2,mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG,g)

Implicit None 
Real(dp), Intent(out) :: g(353) 
Real(dp), Intent(in) :: g1,g2,g3,mHd2,mHu2

Complex(dp), Intent(in) :: Yu(3,3),Yd(3,3),Ye(3,3),Yt(3,3),Ys(3,3),Yz(3,3)            & 
& ,L1,L2,MTM,mue,MZM,MSM,AYu(3,3),AYd(3,3),AYe(3,3),AYt(3,3),AYs(3,3),AYz(3,3)           & 
& ,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2(3,3),ml2(3,3),md2(3,3),mu2(3,3),me2(3,3)              & 
& ,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG

Integer i1, i2, SumI 
 
Iname = Iname +1 
NameOfUnit(Iname) = 'ParametersToG353' 
 
g(1) = g1  
g(2) = g2  
g(3) = g3  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+4) = Real(Yu(i1,i2), dp) 
g(SumI+5) = Aimag(Yu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+22) = Real(Yd(i1,i2), dp) 
g(SumI+23) = Aimag(Yd(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+40) = Real(Ye(i1,i2), dp) 
g(SumI+41) = Aimag(Ye(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+58) = Real(Yt(i1,i2), dp) 
g(SumI+59) = Aimag(Yt(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+76) = Real(Ys(i1,i2), dp) 
g(SumI+77) = Aimag(Ys(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
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
g(118) = Real(mue,dp)  
g(119) = Aimag(mue)  
g(120) = Real(MZM,dp)  
g(121) = Aimag(MZM)  
g(122) = Real(MSM,dp)  
g(123) = Aimag(MSM)  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+124) = Real(AYu(i1,i2), dp) 
g(SumI+125) = Aimag(AYu(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+142) = Real(AYd(i1,i2), dp) 
g(SumI+143) = Aimag(AYd(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+160) = Real(AYe(i1,i2), dp) 
g(SumI+161) = Aimag(AYe(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+178) = Real(AYt(i1,i2), dp) 
g(SumI+179) = Aimag(AYt(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+196) = Real(AYs(i1,i2), dp) 
g(SumI+197) = Aimag(AYs(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+214) = Real(AYz(i1,i2), dp) 
g(SumI+215) = Aimag(AYz(i1,i2)) 
End Do 
End Do 

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
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+244) = Real(mq2(i1,i2), dp) 
g(SumI+245) = Aimag(mq2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+262) = Real(ml2(i1,i2), dp) 
g(SumI+263) = Aimag(ml2(i1,i2)) 
End Do 
End Do 

g(280) = mHd2  
g(281) = mHu2  
Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+282) = Real(md2(i1,i2), dp) 
g(SumI+283) = Aimag(md2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+300) = Real(mu2(i1,i2), dp) 
g(SumI+301) = Aimag(mu2(i1,i2)) 
End Do 
End Do 

Do i1 = 1,3
Do i2 = 1,3
SumI = (i2-1) + (i1-1)*3
SumI = SumI*2 
g(SumI+318) = Real(me2(i1,i2), dp) 
g(SumI+319) = Aimag(me2(i1,i2)) 
End Do 
End Do 

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

 Real(dp) :: g12, g13, g22, g23, g32, g33

Iname = Iname +1 
NameOfUnit(Iname) = 'rge353' 
 
OnlyDiagonal = .Not.GenerationMixing 
q = t 
 
Call GToParameters353(gy,g1,g2,g3,Yu,Yd,Ye,Yt,Ys,Yz,L1,L2,MTM,mue,MZM,MSM,            & 
& AYu,AYd,AYe,AYt,AYs,AYz,AL1,AL2,Amue,AMTM,AMZM,AMSM,mq2,ml2,mHd2,mHu2,md2,             & 
& mu2,me2,mt2,mtb2,ms2,msb2,mz2,mzb2,MassB,MassWB,MassG)

 g12 = g1**2
 g13 = g1*g12
 g22 = g2**2
 g23 = g2*g22
 g32 = g3**2
 g33 = g3*g32

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
 md2Ys = Matmul2(md2,Ys,OnlyDiagonal) 
 md2CYd = Matmul2(md2,Conjg(Yd),OnlyDiagonal) 
 md2CYs = Matmul2(md2,adjYs,OnlyDiagonal) 
 md2CYz = Matmul2(md2,Conjg(Yz),OnlyDiagonal) 
 me2CYe = Matmul2(me2,Conjg(Ye),OnlyDiagonal) 
 ml2Yt = Matmul2(ml2,Yt,OnlyDiagonal) 
 ml2adjYe = Matmul2(ml2,adjYe,OnlyDiagonal) 
 ml2adjYz = Matmul2(ml2,adjYz,OnlyDiagonal) 
 ml2CYt = Matmul2(ml2,adjYt,OnlyDiagonal) 
 mq2adjYd = Matmul2(mq2,adjYd,OnlyDiagonal) 
 mq2adjYu = Matmul2(mq2,adjYu,OnlyDiagonal) 
 mu2CYu = Matmul2(mu2,Conjg(Yu),OnlyDiagonal) 
 YdadjYd = Matmul2(Yd,adjYd,OnlyDiagonal) 
 YeadjYe = Matmul2(Ye,adjYe,OnlyDiagonal) 
 Ysmd2 = Matmul2(Ys,md2,OnlyDiagonal) 
 YsCYs = Matmul2(Ys,adjYs,OnlyDiagonal) 
 Ytml2 = Matmul2(Yt,ml2,OnlyDiagonal) 
 YtCYt = Matmul2(Yt,adjYt,OnlyDiagonal) 
 YuadjYu = Matmul2(Yu,adjYu,OnlyDiagonal) 
 YzadjYz = Matmul2(Yz,adjYz,OnlyDiagonal) 
 AYdadjYd = Matmul2(AYd,adjYd,OnlyDiagonal) 
 AYdadjAYd = Matmul2(AYd,adjAYd,OnlyDiagonal) 
 AYeadjYe = Matmul2(AYe,adjYe,OnlyDiagonal) 
 AYeadjAYe = Matmul2(AYe,adjAYe,OnlyDiagonal) 
 AYsCAYs = Matmul2(AYs,adjAYs,OnlyDiagonal) 
 AYtCAYt = Matmul2(AYt,adjAYt,OnlyDiagonal) 
 AYuadjYu = Matmul2(AYu,adjYu,OnlyDiagonal) 
 AYuadjAYu = Matmul2(AYu,adjAYu,OnlyDiagonal) 
 AYzadjYz = Matmul2(AYz,adjYz,OnlyDiagonal) 
 AYzadjAYz = Matmul2(AYz,adjAYz,OnlyDiagonal) 
 adjYdmd2 = Matmul2(adjYd,md2,OnlyDiagonal) 
 adjYdYd = Matmul2(adjYd,Yd,OnlyDiagonal) 
 adjYdYs = Matmul2(adjYd,Ys,OnlyDiagonal) 
 adjYdYz = Matmul2(adjYd,Yz,OnlyDiagonal) 
 adjYdAYd = Matmul2(adjYd,AYd,OnlyDiagonal) 
 adjYdAYs = Matmul2(adjYd,AYs,OnlyDiagonal) 
 adjYdAYz = Matmul2(adjYd,AYz,OnlyDiagonal) 
 adjYeme2 = Matmul2(adjYe,me2,OnlyDiagonal) 
 adjYeYe = Matmul2(adjYe,Ye,OnlyDiagonal) 
 adjYeAYe = Matmul2(adjYe,AYe,OnlyDiagonal) 
 adjYumu2 = Matmul2(adjYu,mu2,OnlyDiagonal) 
 adjYuYu = Matmul2(adjYu,Yu,OnlyDiagonal) 
 adjYuAYu = Matmul2(adjYu,AYu,OnlyDiagonal) 
 adjYzmd2 = Matmul2(adjYz,md2,OnlyDiagonal) 
 adjYzYd = Matmul2(adjYz,Yd,OnlyDiagonal) 
 adjYzYs = Matmul2(adjYz,Ys,OnlyDiagonal) 
 adjYzYz = Matmul2(adjYz,Yz,OnlyDiagonal) 
 adjYzAYd = Matmul2(adjYz,AYd,OnlyDiagonal) 
 adjYzAYs = Matmul2(adjYz,AYs,OnlyDiagonal) 
 adjYzAYz = Matmul2(adjYz,AYz,OnlyDiagonal) 
 CYdmq2 = Matmul2(Conjg(Yd),mq2,OnlyDiagonal) 
 CYdTYd = Matmul2(Conjg(Yd),Transpose(Yd),OnlyDiagonal) 
 CYdTAYd = Matmul2(Conjg(Yd),Transpose(AYd),OnlyDiagonal) 
 CYeml2 = Matmul2(Conjg(Ye),ml2,OnlyDiagonal) 
 CYeYt = Matmul2(Conjg(Ye),Yt,OnlyDiagonal) 
 CYeAYt = Matmul2(Conjg(Ye),AYt,OnlyDiagonal) 
 CYsmd2 = Matmul2(adjYs,md2,OnlyDiagonal) 
 CYsYd = Matmul2(adjYs,Yd,OnlyDiagonal) 
 CYsYs = Matmul2(adjYs,Ys,OnlyDiagonal) 
 CYsYz = Matmul2(adjYs,Yz,OnlyDiagonal) 
 CYsAYd = Matmul2(adjYs,AYd,OnlyDiagonal) 
 CYsAYs = Matmul2(adjYs,AYs,OnlyDiagonal) 
 CYsAYz = Matmul2(adjYs,AYz,OnlyDiagonal) 
 CYtml2 = Matmul2(adjYt,ml2,OnlyDiagonal) 
 CYtYt = Matmul2(adjYt,Yt,OnlyDiagonal) 
 CYtAYt = Matmul2(adjYt,AYt,OnlyDiagonal) 
 CYumq2 = Matmul2(Conjg(Yu),mq2,OnlyDiagonal) 
 CYzml2 = Matmul2(Conjg(Yz),ml2,OnlyDiagonal) 
 CYzYt = Matmul2(Conjg(Yz),Yt,OnlyDiagonal) 
 CYzAYt = Matmul2(Conjg(Yz),AYt,OnlyDiagonal) 
 CYzTYz = Matmul2(Conjg(Yz),Transpose(Yz),OnlyDiagonal) 
 CYzTAYz = Matmul2(Conjg(Yz),Transpose(AYz),OnlyDiagonal) 
 CAYsAYs = Matmul2(adjAYs,AYs,OnlyDiagonal) 
 CAYtAYt = Matmul2(adjAYt,AYt,OnlyDiagonal) 
 TYdCYd = Matmul2(Transpose(Yd),Conjg(Yd),OnlyDiagonal) 
 TYeCYe = Matmul2(Transpose(Ye),Conjg(Ye),OnlyDiagonal) 
 TYuCYu = Matmul2(Transpose(Yu),Conjg(Yu),OnlyDiagonal) 
 TYzCYz = Matmul2(Transpose(Yz),Conjg(Yz),OnlyDiagonal) 
 TAYdCAYd = Matmul2(Transpose(AYd),Conjg(AYd),OnlyDiagonal) 
 TAYeCAYe = Matmul2(Transpose(AYe),Conjg(AYe),OnlyDiagonal) 
 TAYuCAYu = Matmul2(Transpose(AYu),Conjg(AYu),OnlyDiagonal) 
 TAYzCAYz = Matmul2(Transpose(AYz),Conjg(AYz),OnlyDiagonal) 
 md2YdadjYd = Matmul2(md2,YdadjYd,OnlyDiagonal) 
 md2YsCYs = Matmul2(md2,YsCYs,OnlyDiagonal) 
 md2YzadjYz = Matmul2(md2,YzadjYz,OnlyDiagonal) 
 me2YeadjYe = Matmul2(me2,YeadjYe,OnlyDiagonal) 
 ml2YtCYt = Matmul2(ml2,YtCYt,OnlyDiagonal) 
 ml2TYeCYe = Matmul2(ml2,TYeCYe,OnlyDiagonal) 
 ml2TYzCYz = Matmul2(ml2,TYzCYz,OnlyDiagonal) 
 mq2TYdCYd = Matmul2(mq2,TYdCYd,OnlyDiagonal) 
 mq2TYuCYu = Matmul2(mq2,TYuCYu,OnlyDiagonal) 
 mu2YuadjYu = Matmul2(mu2,YuadjYu,OnlyDiagonal) 
 Ydmq2adjYd = Matmul2(Yd,mq2adjYd,OnlyDiagonal) 
 YdadjYdmd2 = Matmul2(Yd,adjYdmd2,OnlyDiagonal) 
 YdadjYdYd = Matmul2(Yd,adjYdYd,OnlyDiagonal) 
 YdadjYdYs = Matmul2(Yd,adjYdYs,OnlyDiagonal) 
 YdadjYdYz = Matmul2(Yd,adjYdYz,OnlyDiagonal) 
 YdadjYdAYd = Matmul2(Yd,adjYdAYd,OnlyDiagonal) 
 YdadjYdAYs = Matmul2(Yd,adjYdAYs,OnlyDiagonal) 
 YdadjYdAYz = Matmul2(Yd,adjYdAYz,OnlyDiagonal) 
 YdadjYuYu = Matmul2(Yd,adjYuYu,OnlyDiagonal) 
 YdadjYuAYu = Matmul2(Yd,adjYuAYu,OnlyDiagonal) 
 Yeml2adjYe = Matmul2(Ye,ml2adjYe,OnlyDiagonal) 
 YeadjYeme2 = Matmul2(Ye,adjYeme2,OnlyDiagonal) 
 YeadjYeYe = Matmul2(Ye,adjYeYe,OnlyDiagonal) 
 YeadjYeAYe = Matmul2(Ye,adjYeAYe,OnlyDiagonal) 
 YeadjYzYz = Matmul2(Ye,adjYzYz,OnlyDiagonal) 
 YeadjYzAYz = Matmul2(Ye,adjYzAYz,OnlyDiagonal) 
 YeCYtYt = Matmul2(Ye,CYtYt,OnlyDiagonal) 
 YeCYtAYt = Matmul2(Ye,CYtAYt,OnlyDiagonal) 
 Ysmd2CYs = Matmul2(Ys,md2CYs,OnlyDiagonal) 
 YsCYdTYd = Matmul2(Ys,CYdTYd,OnlyDiagonal) 
 YsCYdTAYd = Matmul2(Ys,CYdTAYd,OnlyDiagonal) 
 YsCYsmd2 = Matmul2(Ys,CYsmd2,OnlyDiagonal) 
 YsCYsYd = Matmul2(Ys,CYsYd,OnlyDiagonal) 
 YsCYsYs = Matmul2(Ys,CYsYs,OnlyDiagonal) 
 YsCYsYz = Matmul2(Ys,CYsYz,OnlyDiagonal) 
 YsCYsAYd = Matmul2(Ys,CYsAYd,OnlyDiagonal) 
 YsCYsAYs = Matmul2(Ys,CYsAYs,OnlyDiagonal) 
 YsCYsAYz = Matmul2(Ys,CYsAYz,OnlyDiagonal) 
 YsCYzTYz = Matmul2(Ys,CYzTYz,OnlyDiagonal) 
 YsCYzTAYz = Matmul2(Ys,CYzTAYz,OnlyDiagonal) 
 Ytml2CYt = Matmul2(Yt,ml2CYt,OnlyDiagonal) 
 YtadjYeYe = Matmul2(Yt,adjYeYe,OnlyDiagonal) 
 YtadjYeAYe = Matmul2(Yt,adjYeAYe,OnlyDiagonal) 
 YtadjYzYz = Matmul2(Yt,adjYzYz,OnlyDiagonal) 
 YtadjYzAYz = Matmul2(Yt,adjYzAYz,OnlyDiagonal) 
 YtCYtml2 = Matmul2(Yt,CYtml2,OnlyDiagonal) 
 YtCYtYt = Matmul2(Yt,CYtYt,OnlyDiagonal) 
 YtCYtAYt = Matmul2(Yt,CYtAYt,OnlyDiagonal) 
 Yumq2adjYu = Matmul2(Yu,mq2adjYu,OnlyDiagonal) 
 YuadjYdYd = Matmul2(Yu,adjYdYd,OnlyDiagonal) 
 YuadjYdAYd = Matmul2(Yu,adjYdAYd,OnlyDiagonal) 
 YuadjYumu2 = Matmul2(Yu,adjYumu2,OnlyDiagonal) 
 YuadjYuYu = Matmul2(Yu,adjYuYu,OnlyDiagonal) 
 YuadjYuAYu = Matmul2(Yu,adjYuAYu,OnlyDiagonal) 
 Yzml2adjYz = Matmul2(Yz,ml2adjYz,OnlyDiagonal) 
 YzadjYeYe = Matmul2(Yz,adjYeYe,OnlyDiagonal) 
 YzadjYeAYe = Matmul2(Yz,adjYeAYe,OnlyDiagonal) 
 YzadjYzmd2 = Matmul2(Yz,adjYzmd2,OnlyDiagonal) 
 YzadjYzYd = Matmul2(Yz,adjYzYd,OnlyDiagonal) 
 YzadjYzYs = Matmul2(Yz,adjYzYs,OnlyDiagonal) 
 YzadjYzYz = Matmul2(Yz,adjYzYz,OnlyDiagonal) 
 YzadjYzAYd = Matmul2(Yz,adjYzAYd,OnlyDiagonal) 
 YzadjYzAYs = Matmul2(Yz,adjYzAYs,OnlyDiagonal) 
 YzadjYzAYz = Matmul2(Yz,adjYzAYz,OnlyDiagonal) 
 YzCYtYt = Matmul2(Yz,CYtYt,OnlyDiagonal) 
 YzCYtAYt = Matmul2(Yz,CYtAYt,OnlyDiagonal) 
 AYdadjYdYd = Matmul2(AYd,adjYdYd,OnlyDiagonal) 
 AYdadjYdYs = Matmul2(AYd,adjYdYs,OnlyDiagonal) 
 AYdadjYdYz = Matmul2(AYd,adjYdYz,OnlyDiagonal) 
 AYdadjYuYu = Matmul2(AYd,adjYuYu,OnlyDiagonal) 
 AYeadjYeYe = Matmul2(AYe,adjYeYe,OnlyDiagonal) 
 AYeadjYzYz = Matmul2(AYe,adjYzYz,OnlyDiagonal) 
 AYeCYtYt = Matmul2(AYe,CYtYt,OnlyDiagonal) 
 AYsCYdTYd = Matmul2(AYs,CYdTYd,OnlyDiagonal) 
 AYsCYsYd = Matmul2(AYs,CYsYd,OnlyDiagonal) 
 AYsCYsYs = Matmul2(AYs,CYsYs,OnlyDiagonal) 
 AYsCYsYz = Matmul2(AYs,CYsYz,OnlyDiagonal) 
 AYsCYzTYz = Matmul2(AYs,CYzTYz,OnlyDiagonal) 
 AYtadjYeYe = Matmul2(AYt,adjYeYe,OnlyDiagonal) 
 AYtadjYzYz = Matmul2(AYt,adjYzYz,OnlyDiagonal) 
 AYtCYtYt = Matmul2(AYt,CYtYt,OnlyDiagonal) 
 AYuadjYdYd = Matmul2(AYu,adjYdYd,OnlyDiagonal) 
 AYuadjYuYu = Matmul2(AYu,adjYuYu,OnlyDiagonal) 
 AYzadjYeYe = Matmul2(AYz,adjYeYe,OnlyDiagonal) 
 AYzadjYzYd = Matmul2(AYz,adjYzYd,OnlyDiagonal) 
 AYzadjYzYs = Matmul2(AYz,adjYzYs,OnlyDiagonal) 
 AYzadjYzYz = Matmul2(AYz,adjYzYz,OnlyDiagonal) 
 AYzCYtYt = Matmul2(AYz,CYtYt,OnlyDiagonal) 
 CYsmd2Ys = Matmul2(adjYs,md2Ys,OnlyDiagonal) 
 CYsYsmd2 = Matmul2(adjYs,Ysmd2,OnlyDiagonal) 
 CYtml2Yt = Matmul2(adjYt,ml2Yt,OnlyDiagonal) 
 CYtYtml2 = Matmul2(adjYt,Ytml2,OnlyDiagonal) 
 TYdmd2CYd = Matmul2(Transpose(Yd),md2CYd,OnlyDiagonal) 
 TYdCYdmq2 = Matmul2(Transpose(Yd),CYdmq2,OnlyDiagonal) 
 TYeme2CYe = Matmul2(Transpose(Ye),me2CYe,OnlyDiagonal) 
 TYeCYeml2 = Matmul2(Transpose(Ye),CYeml2,OnlyDiagonal) 
 TYeCYeYt = Matmul2(Transpose(Ye),CYeYt,OnlyDiagonal) 
 TYeCYeAYt = Matmul2(Transpose(Ye),CYeAYt,OnlyDiagonal) 
 TYumu2CYu = Matmul2(Transpose(Yu),mu2CYu,OnlyDiagonal) 
 TYuCYumq2 = Matmul2(Transpose(Yu),CYumq2,OnlyDiagonal) 
 TYzmd2CYz = Matmul2(Transpose(Yz),md2CYz,OnlyDiagonal) 
 TYzCYzml2 = Matmul2(Transpose(Yz),CYzml2,OnlyDiagonal) 
 TYzCYzYt = Matmul2(Transpose(Yz),CYzYt,OnlyDiagonal) 
 TYzCYzAYt = Matmul2(Transpose(Yz),CYzAYt,OnlyDiagonal) 
 TAYeCYeYt = Matmul2(Transpose(AYe),CYeYt,OnlyDiagonal) 
 TAYzCYzYt = Matmul2(Transpose(AYz),CYzYt,OnlyDiagonal) 
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
 md2Yz = Matmul2(md2,Yz,OnlyDiagonal) 
 me2Ye = Matmul2(me2,Ye,OnlyDiagonal) 
 YdadjYu = Matmul2(Yd,adjYu,OnlyDiagonal) 
 YdadjAYd = Matmul2(Yd,adjAYd,OnlyDiagonal) 
 YdadjAYu = Matmul2(Yd,adjAYu,OnlyDiagonal) 
 Yeml2 = Matmul2(Ye,ml2,OnlyDiagonal) 
 YeadjYz = Matmul2(Ye,adjYz,OnlyDiagonal) 
 YeadjAYe = Matmul2(Ye,adjAYe,OnlyDiagonal) 
 YeadjAYz = Matmul2(Ye,adjAYz,OnlyDiagonal) 
 YeCYt = Matmul2(Ye,adjYt,OnlyDiagonal) 
 YeCAYt = Matmul2(Ye,adjAYt,OnlyDiagonal) 
 YsCYd = Matmul2(Ys,Conjg(Yd),OnlyDiagonal) 
 YsCYz = Matmul2(Ys,Conjg(Yz),OnlyDiagonal) 
 YsCAYd = Matmul2(Ys,Conjg(AYd),OnlyDiagonal) 
 YsCAYs = Matmul2(Ys,adjAYs,OnlyDiagonal) 
 YsCAYz = Matmul2(Ys,Conjg(AYz),OnlyDiagonal) 
 YtadjYe = Matmul2(Yt,adjYe,OnlyDiagonal) 
 YtadjYz = Matmul2(Yt,adjYz,OnlyDiagonal) 
 YtadjAYe = Matmul2(Yt,adjAYe,OnlyDiagonal) 
 YtadjAYz = Matmul2(Yt,adjAYz,OnlyDiagonal) 
 YtCAYt = Matmul2(Yt,adjAYt,OnlyDiagonal) 
 YuadjYd = Matmul2(Yu,adjYd,OnlyDiagonal) 
 YuadjAYd = Matmul2(Yu,adjAYd,OnlyDiagonal) 
 YuadjAYu = Matmul2(Yu,adjAYu,OnlyDiagonal) 
 Yzml2 = Matmul2(Yz,ml2,OnlyDiagonal) 
 YzadjYe = Matmul2(Yz,adjYe,OnlyDiagonal) 
 YzadjAYe = Matmul2(Yz,adjAYe,OnlyDiagonal) 
 YzadjAYz = Matmul2(Yz,adjAYz,OnlyDiagonal) 
 YzCYt = Matmul2(Yz,adjYt,OnlyDiagonal) 
 YzCAYt = Matmul2(Yz,adjAYt,OnlyDiagonal) 
 AYdadjYu = Matmul2(AYd,adjYu,OnlyDiagonal) 
 AYdadjAYu = Matmul2(AYd,adjAYu,OnlyDiagonal) 
 AYeadjYz = Matmul2(AYe,adjYz,OnlyDiagonal) 
 AYeadjAYz = Matmul2(AYe,adjAYz,OnlyDiagonal) 
 AYeCYt = Matmul2(AYe,adjYt,OnlyDiagonal) 
 AYeCAYt = Matmul2(AYe,adjAYt,OnlyDiagonal) 
 AYsCYd = Matmul2(AYs,Conjg(Yd),OnlyDiagonal) 
 AYsCYs = Matmul2(AYs,adjYs,OnlyDiagonal) 
 AYsCYz = Matmul2(AYs,Conjg(Yz),OnlyDiagonal) 
 AYsCAYd = Matmul2(AYs,Conjg(AYd),OnlyDiagonal) 
 AYsCAYz = Matmul2(AYs,Conjg(AYz),OnlyDiagonal) 
 AYtadjYe = Matmul2(AYt,adjYe,OnlyDiagonal) 
 AYtadjYz = Matmul2(AYt,adjYz,OnlyDiagonal) 
 AYtadjAYe = Matmul2(AYt,adjAYe,OnlyDiagonal) 
 AYtadjAYz = Matmul2(AYt,adjAYz,OnlyDiagonal) 
 AYtCYt = Matmul2(AYt,adjYt,OnlyDiagonal) 
 AYuadjYd = Matmul2(AYu,adjYd,OnlyDiagonal) 
 AYuadjAYd = Matmul2(AYu,adjAYd,OnlyDiagonal) 
 AYzadjYe = Matmul2(AYz,adjYe,OnlyDiagonal) 
 AYzadjAYe = Matmul2(AYz,adjAYe,OnlyDiagonal) 
 AYzCYt = Matmul2(AYz,adjYt,OnlyDiagonal) 
 AYzCAYt = Matmul2(AYz,adjAYt,OnlyDiagonal) 
 adjAYdYs = Matmul2(adjAYd,Ys,OnlyDiagonal) 
 adjAYdAYs = Matmul2(adjAYd,AYs,OnlyDiagonal) 
 adjAYeYe = Matmul2(adjAYe,Ye,OnlyDiagonal) 
 adjAYeAYe = Matmul2(adjAYe,AYe,OnlyDiagonal) 
 adjAYzYs = Matmul2(adjAYz,Ys,OnlyDiagonal) 
 adjAYzYz = Matmul2(adjAYz,Yz,OnlyDiagonal) 
 adjAYzAYs = Matmul2(adjAYz,AYs,OnlyDiagonal) 
 adjAYzAYz = Matmul2(adjAYz,AYz,OnlyDiagonal) 
 CYeTYz = Matmul2(Conjg(Ye),Transpose(Yz),OnlyDiagonal) 
 CYeTAYz = Matmul2(Conjg(Ye),Transpose(AYz),OnlyDiagonal) 
 CYtTYz = Matmul2(adjYt,Transpose(Yz),OnlyDiagonal) 
 CYtTAYz = Matmul2(adjYt,Transpose(AYz),OnlyDiagonal) 
 CYuTYd = Matmul2(Conjg(Yu),Transpose(Yd),OnlyDiagonal) 
 CYuTAYd = Matmul2(Conjg(Yu),Transpose(AYd),OnlyDiagonal) 
 CAYsYs = Matmul2(adjAYs,Ys,OnlyDiagonal) 
 CAYtYt = Matmul2(adjAYt,Yt,OnlyDiagonal) 
 TYdCYs = Matmul2(Transpose(Yd),adjYs,OnlyDiagonal) 
 TYdCYz = Matmul2(Transpose(Yd),Conjg(Yz),OnlyDiagonal) 
 TYdCAYd = Matmul2(Transpose(Yd),Conjg(AYd),OnlyDiagonal) 
 TYdCAYs = Matmul2(Transpose(Yd),adjAYs,OnlyDiagonal) 
 TYdCAYz = Matmul2(Transpose(Yd),Conjg(AYz),OnlyDiagonal) 
 TYeCAYe = Matmul2(Transpose(Ye),Conjg(AYe),OnlyDiagonal) 
 TYuCAYu = Matmul2(Transpose(Yu),Conjg(AYu),OnlyDiagonal) 
 TYzCYd = Matmul2(Transpose(Yz),Conjg(Yd),OnlyDiagonal) 
 TYzCYs = Matmul2(Transpose(Yz),adjYs,OnlyDiagonal) 
 TYzCAYd = Matmul2(Transpose(Yz),Conjg(AYd),OnlyDiagonal) 
 TYzCAYs = Matmul2(Transpose(Yz),adjAYs,OnlyDiagonal) 
 TYzCAYz = Matmul2(Transpose(Yz),Conjg(AYz),OnlyDiagonal) 
 TAYdCYd = Matmul2(Transpose(AYd),Conjg(Yd),OnlyDiagonal) 
 TAYdCYs = Matmul2(Transpose(AYd),adjYs,OnlyDiagonal) 
 TAYdCYz = Matmul2(Transpose(AYd),Conjg(Yz),OnlyDiagonal) 
 TAYdCAYs = Matmul2(Transpose(AYd),adjAYs,OnlyDiagonal) 
 TAYdCAYz = Matmul2(Transpose(AYd),Conjg(AYz),OnlyDiagonal) 
 TAYeCYe = Matmul2(Transpose(AYe),Conjg(Ye),OnlyDiagonal) 
 TAYuCYu = Matmul2(Transpose(AYu),Conjg(Yu),OnlyDiagonal) 
 TAYzCYd = Matmul2(Transpose(AYz),Conjg(Yd),OnlyDiagonal) 
 TAYzCYs = Matmul2(Transpose(AYz),adjYs,OnlyDiagonal) 
 TAYzCYz = Matmul2(Transpose(AYz),Conjg(Yz),OnlyDiagonal) 
 TAYzCAYd = Matmul2(Transpose(AYz),Conjg(AYd),OnlyDiagonal) 
 TAYzCAYs = Matmul2(Transpose(AYz),adjAYs,OnlyDiagonal) 
 md2YdadjYu = Matmul2(md2,YdadjYu,OnlyDiagonal) 
 md2YsCYd = Matmul2(md2,YsCYd,OnlyDiagonal) 
 md2YsCYz = Matmul2(md2,YsCYz,OnlyDiagonal) 
 md2YzadjYe = Matmul2(md2,YzadjYe,OnlyDiagonal) 
 md2YzCYt = Matmul2(md2,YzCYt,OnlyDiagonal) 
 md2CYsYs = Matmul2(md2,CYsYs,OnlyDiagonal) 
 me2YeadjYz = Matmul2(me2,YeadjYz,OnlyDiagonal) 
 me2YeCYt = Matmul2(me2,YeCYt,OnlyDiagonal) 
 ml2YtadjYe = Matmul2(ml2,YtadjYe,OnlyDiagonal) 
 ml2YtadjYz = Matmul2(ml2,YtadjYz,OnlyDiagonal) 
 ml2adjYeYe = Matmul2(ml2,adjYeYe,OnlyDiagonal) 
 ml2adjYzYs = Matmul2(ml2,adjYzYs,OnlyDiagonal) 
 ml2adjYzYz = Matmul2(ml2,adjYzYz,OnlyDiagonal) 
 ml2CYtYt = Matmul2(ml2,CYtYt,OnlyDiagonal) 
 ml2TYzCYd = Matmul2(ml2,TYzCYd,OnlyDiagonal) 
 ml2TYzCYs = Matmul2(ml2,TYzCYs,OnlyDiagonal) 
 mq2adjYdYs = Matmul2(mq2,adjYdYs,OnlyDiagonal) 
 mq2TYdCYs = Matmul2(mq2,TYdCYs,OnlyDiagonal) 
 mq2TYdCYz = Matmul2(mq2,TYdCYz,OnlyDiagonal) 
 mu2YuadjYd = Matmul2(mu2,YuadjYd,OnlyDiagonal) 
 Ydmq2adjYu = Matmul2(Yd,mq2adjYu,OnlyDiagonal) 
 YdadjYumu2 = Matmul2(Yd,adjYumu2,OnlyDiagonal) 
 YdadjAYdAYs = Matmul2(Yd,adjAYdAYs,OnlyDiagonal) 
 Yeml2adjYz = Matmul2(Ye,ml2adjYz,OnlyDiagonal) 
 Yeml2CYt = Matmul2(Ye,ml2CYt,OnlyDiagonal) 
 YeadjYzmd2 = Matmul2(Ye,adjYzmd2,OnlyDiagonal) 
 YeadjYzYd = Matmul2(Ye,adjYzYd,OnlyDiagonal) 
 YeadjYzYs = Matmul2(Ye,adjYzYs,OnlyDiagonal) 
 YeadjYzAYd = Matmul2(Ye,adjYzAYd,OnlyDiagonal) 
 YeadjYzAYs = Matmul2(Ye,adjYzAYs,OnlyDiagonal) 
 YeCYtml2 = Matmul2(Ye,CYtml2,OnlyDiagonal) 
 Ysmd2CYd = Matmul2(Ys,md2CYd,OnlyDiagonal) 
 Ysmd2CYz = Matmul2(Ys,md2CYz,OnlyDiagonal) 
 YsCYdmq2 = Matmul2(Ys,CYdmq2,OnlyDiagonal) 
 YsCYzml2 = Matmul2(Ys,CYzml2,OnlyDiagonal) 
 YsCYzYt = Matmul2(Ys,CYzYt,OnlyDiagonal) 
 YsCYzAYt = Matmul2(Ys,CYzAYt,OnlyDiagonal) 
 YsCAYsAYs = Matmul2(Ys,CAYsAYs,OnlyDiagonal) 
 Ytml2adjYe = Matmul2(Yt,ml2adjYe,OnlyDiagonal) 
 Ytml2adjYz = Matmul2(Yt,ml2adjYz,OnlyDiagonal) 
 YtadjYeme2 = Matmul2(Yt,adjYeme2,OnlyDiagonal) 
 YtadjYzmd2 = Matmul2(Yt,adjYzmd2,OnlyDiagonal) 
 YtadjYzYd = Matmul2(Yt,adjYzYd,OnlyDiagonal) 
 YtadjYzYs = Matmul2(Yt,adjYzYs,OnlyDiagonal) 
 YtadjYzAYd = Matmul2(Yt,adjYzAYd,OnlyDiagonal) 
 YtadjYzAYs = Matmul2(Yt,adjYzAYs,OnlyDiagonal) 
 YtadjAYeAYe = Matmul2(Yt,adjAYeAYe,OnlyDiagonal) 
 YtadjAYzAYz = Matmul2(Yt,adjAYzAYz,OnlyDiagonal) 
 YtCYtTYz = Matmul2(Yt,CYtTYz,OnlyDiagonal) 
 YtCYtTAYz = Matmul2(Yt,CYtTAYz,OnlyDiagonal) 
 YtCAYtAYt = Matmul2(Yt,CAYtAYt,OnlyDiagonal) 
 Yumq2adjYd = Matmul2(Yu,mq2adjYd,OnlyDiagonal) 
 YuadjYdmd2 = Matmul2(Yu,adjYdmd2,OnlyDiagonal) 
 YuadjYdYs = Matmul2(Yu,adjYdYs,OnlyDiagonal) 
 YuadjYdYz = Matmul2(Yu,adjYdYz,OnlyDiagonal) 
 YuadjYdAYs = Matmul2(Yu,adjYdAYs,OnlyDiagonal) 
 YuadjYdAYz = Matmul2(Yu,adjYdAYz,OnlyDiagonal) 
 Yzml2adjYe = Matmul2(Yz,ml2adjYe,OnlyDiagonal) 
 Yzml2CYt = Matmul2(Yz,ml2CYt,OnlyDiagonal) 
 YzadjYeme2 = Matmul2(Yz,adjYeme2,OnlyDiagonal) 
 YzadjAYzAYs = Matmul2(Yz,adjAYzAYs,OnlyDiagonal) 
 YzCYtml2 = Matmul2(Yz,CYtml2,OnlyDiagonal) 
 AYdadjAYdYs = Matmul2(AYd,adjAYdYs,OnlyDiagonal) 
 AYeadjYzYd = Matmul2(AYe,adjYzYd,OnlyDiagonal) 
 AYeadjYzYs = Matmul2(AYe,adjYzYs,OnlyDiagonal) 
 AYsCYzYt = Matmul2(AYs,CYzYt,OnlyDiagonal) 
 AYsCAYsYs = Matmul2(AYs,CAYsYs,OnlyDiagonal) 
 AYtadjYzYd = Matmul2(AYt,adjYzYd,OnlyDiagonal) 
 AYtadjYzYs = Matmul2(AYt,adjYzYs,OnlyDiagonal) 
 AYtadjAYeYe = Matmul2(AYt,adjAYeYe,OnlyDiagonal) 
 AYtadjAYzYz = Matmul2(AYt,adjAYzYz,OnlyDiagonal) 
 AYtCYtTYz = Matmul2(AYt,CYtTYz,OnlyDiagonal) 
 AYtCAYtYt = Matmul2(AYt,CAYtYt,OnlyDiagonal) 
 AYuadjYdYs = Matmul2(AYu,adjYdYs,OnlyDiagonal) 
 AYuadjYdYz = Matmul2(AYu,adjYdYz,OnlyDiagonal) 
 AYzadjAYzYs = Matmul2(AYz,adjAYzYs,OnlyDiagonal) 
 adjYdmd2Ys = Matmul2(adjYd,md2Ys,OnlyDiagonal) 
 adjYdYdadjYd = Matmul2(adjYd,YdadjYd,OnlyDiagonal) 
 adjYdYdadjYu = Matmul2(adjYd,YdadjYu,OnlyDiagonal) 
 adjYdYdadjAYd = Matmul2(adjYd,YdadjAYd,OnlyDiagonal) 
 adjYdYdadjAYu = Matmul2(adjYd,YdadjAYu,OnlyDiagonal) 
 adjYdYsmd2 = Matmul2(adjYd,Ysmd2,OnlyDiagonal) 
 adjYdYsCYs = Matmul2(adjYd,YsCYs,OnlyDiagonal) 
 adjYdYzadjYz = Matmul2(adjYd,YzadjYz,OnlyDiagonal) 
 adjYdAYdadjYd = Matmul2(adjYd,AYdadjYd,OnlyDiagonal) 
 adjYdAYdadjYu = Matmul2(adjYd,AYdadjYu,OnlyDiagonal) 
 adjYdAYdadjAYd = Matmul2(adjYd,AYdadjAYd,OnlyDiagonal) 
 adjYdAYdadjAYu = Matmul2(adjYd,AYdadjAYu,OnlyDiagonal) 
 adjYdAYsCYs = Matmul2(adjYd,AYsCYs,OnlyDiagonal) 
 adjYdAYsCAYs = Matmul2(adjYd,AYsCAYs,OnlyDiagonal) 
 adjYdAYzadjYz = Matmul2(adjYd,AYzadjYz,OnlyDiagonal) 
 adjYdAYzadjAYz = Matmul2(adjYd,AYzadjAYz,OnlyDiagonal) 
 adjYeme2Ye = Matmul2(adjYe,me2Ye,OnlyDiagonal) 
 adjYeYeml2 = Matmul2(adjYe,Yeml2,OnlyDiagonal) 
 adjYeYeadjYe = Matmul2(adjYe,YeadjYe,OnlyDiagonal) 
 adjYeYeadjYz = Matmul2(adjYe,YeadjYz,OnlyDiagonal) 
 adjYeYeadjAYe = Matmul2(adjYe,YeadjAYe,OnlyDiagonal) 
 adjYeYeadjAYz = Matmul2(adjYe,YeadjAYz,OnlyDiagonal) 
 adjYeYeCYt = Matmul2(adjYe,YeCYt,OnlyDiagonal) 
 adjYeYeCAYt = Matmul2(adjYe,YeCAYt,OnlyDiagonal) 
 adjYeAYeadjYe = Matmul2(adjYe,AYeadjYe,OnlyDiagonal) 
 adjYeAYeadjYz = Matmul2(adjYe,AYeadjYz,OnlyDiagonal) 
 adjYeAYeadjAYe = Matmul2(adjYe,AYeadjAYe,OnlyDiagonal) 
 adjYeAYeadjAYz = Matmul2(adjYe,AYeadjAYz,OnlyDiagonal) 
 adjYeAYeCYt = Matmul2(adjYe,AYeCYt,OnlyDiagonal) 
 adjYeAYeCAYt = Matmul2(adjYe,AYeCAYt,OnlyDiagonal) 
 adjYuYuadjYd = Matmul2(adjYu,YuadjYd,OnlyDiagonal) 
 adjYuYuadjYu = Matmul2(adjYu,YuadjYu,OnlyDiagonal) 
 adjYuYuadjAYd = Matmul2(adjYu,YuadjAYd,OnlyDiagonal) 
 adjYuYuadjAYu = Matmul2(adjYu,YuadjAYu,OnlyDiagonal) 
 adjYuAYuadjYd = Matmul2(adjYu,AYuadjYd,OnlyDiagonal) 
 adjYuAYuadjYu = Matmul2(adjYu,AYuadjYu,OnlyDiagonal) 
 adjYuAYuadjAYd = Matmul2(adjYu,AYuadjAYd,OnlyDiagonal) 
 adjYuAYuadjAYu = Matmul2(adjYu,AYuadjAYu,OnlyDiagonal) 
 adjYzmd2Ys = Matmul2(adjYz,md2Ys,OnlyDiagonal) 
 adjYzmd2Yz = Matmul2(adjYz,md2Yz,OnlyDiagonal) 
 adjYzYdadjYd = Matmul2(adjYz,YdadjYd,OnlyDiagonal) 
 adjYzYsmd2 = Matmul2(adjYz,Ysmd2,OnlyDiagonal) 
 adjYzYsCYs = Matmul2(adjYz,YsCYs,OnlyDiagonal) 
 adjYzYzml2 = Matmul2(adjYz,Yzml2,OnlyDiagonal) 
 adjYzYzadjYe = Matmul2(adjYz,YzadjYe,OnlyDiagonal) 
 adjYzYzadjYz = Matmul2(adjYz,YzadjYz,OnlyDiagonal) 
 adjYzYzadjAYe = Matmul2(adjYz,YzadjAYe,OnlyDiagonal) 
 adjYzYzadjAYz = Matmul2(adjYz,YzadjAYz,OnlyDiagonal) 
 adjYzYzCYt = Matmul2(adjYz,YzCYt,OnlyDiagonal) 
 adjYzYzCAYt = Matmul2(adjYz,YzCAYt,OnlyDiagonal) 
 adjYzAYdadjYd = Matmul2(adjYz,AYdadjYd,OnlyDiagonal) 
 adjYzAYdadjAYd = Matmul2(adjYz,AYdadjAYd,OnlyDiagonal) 
 adjYzAYsCYs = Matmul2(adjYz,AYsCYs,OnlyDiagonal) 
 adjYzAYsCAYs = Matmul2(adjYz,AYsCAYs,OnlyDiagonal) 
 adjYzAYzadjYe = Matmul2(adjYz,AYzadjYe,OnlyDiagonal) 
 adjYzAYzadjYz = Matmul2(adjYz,AYzadjYz,OnlyDiagonal) 
 adjYzAYzadjAYe = Matmul2(adjYz,AYzadjAYe,OnlyDiagonal) 
 adjYzAYzadjAYz = Matmul2(adjYz,AYzadjAYz,OnlyDiagonal) 
 adjYzAYzCYt = Matmul2(adjYz,AYzCYt,OnlyDiagonal) 
 adjYzAYzCAYt = Matmul2(adjYz,AYzCAYt,OnlyDiagonal) 
 adjAYdYdadjYd = Matmul2(adjAYd,YdadjYd,OnlyDiagonal) 
 adjAYdYdadjYu = Matmul2(adjAYd,YdadjYu,OnlyDiagonal) 
 adjAYdAYdadjYd = Matmul2(adjAYd,AYdadjYd,OnlyDiagonal) 
 adjAYdAYdadjYu = Matmul2(adjAYd,AYdadjYu,OnlyDiagonal) 
 adjAYdAYsCYs = Matmul2(adjAYd,AYsCYs,OnlyDiagonal) 
 adjAYdAYzadjYz = Matmul2(adjAYd,AYzadjYz,OnlyDiagonal) 
 adjAYeYeadjYe = Matmul2(adjAYe,YeadjYe,OnlyDiagonal) 
 adjAYeYeadjYz = Matmul2(adjAYe,YeadjYz,OnlyDiagonal) 
 adjAYeYeCYt = Matmul2(adjAYe,YeCYt,OnlyDiagonal) 
 adjAYeAYeadjYe = Matmul2(adjAYe,AYeadjYe,OnlyDiagonal) 
 adjAYeAYeadjYz = Matmul2(adjAYe,AYeadjYz,OnlyDiagonal) 
 adjAYeAYeCYt = Matmul2(adjAYe,AYeCYt,OnlyDiagonal) 
 adjAYuYuadjYd = Matmul2(adjAYu,YuadjYd,OnlyDiagonal) 
 adjAYuYuadjYu = Matmul2(adjAYu,YuadjYu,OnlyDiagonal) 
 adjAYuAYuadjYd = Matmul2(adjAYu,AYuadjYd,OnlyDiagonal) 
 adjAYuAYuadjYu = Matmul2(adjAYu,AYuadjYu,OnlyDiagonal) 
 adjAYzYzadjYe = Matmul2(adjAYz,YzadjYe,OnlyDiagonal) 
 adjAYzYzadjYz = Matmul2(adjAYz,YzadjYz,OnlyDiagonal) 
 adjAYzYzCYt = Matmul2(adjAYz,YzCYt,OnlyDiagonal) 
 adjAYzAYdadjYd = Matmul2(adjAYz,AYdadjYd,OnlyDiagonal) 
 adjAYzAYsCYs = Matmul2(adjAYz,AYsCYs,OnlyDiagonal) 
 adjAYzAYzadjYe = Matmul2(adjAYz,AYzadjYe,OnlyDiagonal) 
 adjAYzAYzadjYz = Matmul2(adjAYz,AYzadjYz,OnlyDiagonal) 
 adjAYzAYzCYt = Matmul2(adjAYz,AYzCYt,OnlyDiagonal) 
 CYdTYdCYd = Matmul2(Conjg(Yd),TYdCYd,OnlyDiagonal) 
 CYdTYdCYs = Matmul2(Conjg(Yd),TYdCYs,OnlyDiagonal) 
 CYdTYdCYz = Matmul2(Conjg(Yd),TYdCYz,OnlyDiagonal) 
 CYdTYdCAYd = Matmul2(Conjg(Yd),TYdCAYd,OnlyDiagonal) 
 CYdTYdCAYs = Matmul2(Conjg(Yd),TYdCAYs,OnlyDiagonal) 
 CYdTYdCAYz = Matmul2(Conjg(Yd),TYdCAYz,OnlyDiagonal) 
 CYdTAYdCAYd = Matmul2(Conjg(Yd),TAYdCAYd,OnlyDiagonal) 
 CYdTAYdCAYs = Matmul2(Conjg(Yd),TAYdCAYs,OnlyDiagonal) 
 CYdTAYdCAYz = Matmul2(Conjg(Yd),TAYdCAYz,OnlyDiagonal) 
 CYeTYeCYe = Matmul2(Conjg(Ye),TYeCYe,OnlyDiagonal) 
 CYeTYeCAYe = Matmul2(Conjg(Ye),TYeCAYe,OnlyDiagonal) 
 CYeTAYeCAYe = Matmul2(Conjg(Ye),TAYeCAYe,OnlyDiagonal) 
 CYsYdadjYd = Matmul2(adjYs,YdadjYd,OnlyDiagonal) 
 CYsYsCYd = Matmul2(adjYs,YsCYd,OnlyDiagonal) 
 CYsYsCYs = Matmul2(adjYs,YsCYs,OnlyDiagonal) 
 CYsYsCYz = Matmul2(adjYs,YsCYz,OnlyDiagonal) 
 CYsYsCAYd = Matmul2(adjYs,YsCAYd,OnlyDiagonal) 
 CYsYsCAYs = Matmul2(adjYs,YsCAYs,OnlyDiagonal) 
 CYsYsCAYz = Matmul2(adjYs,YsCAYz,OnlyDiagonal) 
 CYsYzadjYz = Matmul2(adjYs,YzadjYz,OnlyDiagonal) 
 CYsAYdadjYd = Matmul2(adjYs,AYdadjYd,OnlyDiagonal) 
 CYsAYdadjAYd = Matmul2(adjYs,AYdadjAYd,OnlyDiagonal) 
 CYsAYsCYs = Matmul2(adjYs,AYsCYs,OnlyDiagonal) 
 CYsAYsCAYd = Matmul2(adjYs,AYsCAYd,OnlyDiagonal) 
 CYsAYsCAYs = Matmul2(adjYs,AYsCAYs,OnlyDiagonal) 
 CYsAYsCAYz = Matmul2(adjYs,AYsCAYz,OnlyDiagonal) 
 CYsAYzadjYz = Matmul2(adjYs,AYzadjYz,OnlyDiagonal) 
 CYsAYzadjAYz = Matmul2(adjYs,AYzadjAYz,OnlyDiagonal) 
 CYtYtadjYe = Matmul2(adjYt,YtadjYe,OnlyDiagonal) 
 CYtYtadjYz = Matmul2(adjYt,YtadjYz,OnlyDiagonal) 
 CYtYtadjAYe = Matmul2(adjYt,YtadjAYe,OnlyDiagonal) 
 CYtYtadjAYz = Matmul2(adjYt,YtadjAYz,OnlyDiagonal) 
 CYtYtCYt = Matmul2(adjYt,YtCYt,OnlyDiagonal) 
 CYtYtCAYt = Matmul2(adjYt,YtCAYt,OnlyDiagonal) 
 CYtAYtadjYe = Matmul2(adjYt,AYtadjYe,OnlyDiagonal) 
 CYtAYtadjYz = Matmul2(adjYt,AYtadjYz,OnlyDiagonal) 
 CYtAYtadjAYe = Matmul2(adjYt,AYtadjAYe,OnlyDiagonal) 
 CYtAYtadjAYz = Matmul2(adjYt,AYtadjAYz,OnlyDiagonal) 
 CYtAYtCYt = Matmul2(adjYt,AYtCYt,OnlyDiagonal) 
 CYtAYtCAYt = Matmul2(adjYt,AYtCAYt,OnlyDiagonal) 
 CYuTYuCYu = Matmul2(Conjg(Yu),TYuCYu,OnlyDiagonal) 
 CYuTYuCAYu = Matmul2(Conjg(Yu),TYuCAYu,OnlyDiagonal) 
 CYuTAYuCAYu = Matmul2(Conjg(Yu),TAYuCAYu,OnlyDiagonal) 
 CYzTYzCYd = Matmul2(Conjg(Yz),TYzCYd,OnlyDiagonal) 
 CYzTYzCYs = Matmul2(Conjg(Yz),TYzCYs,OnlyDiagonal) 
 CYzTYzCYz = Matmul2(Conjg(Yz),TYzCYz,OnlyDiagonal) 
 CYzTYzCAYd = Matmul2(Conjg(Yz),TYzCAYd,OnlyDiagonal) 
 CYzTYzCAYs = Matmul2(Conjg(Yz),TYzCAYs,OnlyDiagonal) 
 CYzTYzCAYz = Matmul2(Conjg(Yz),TYzCAYz,OnlyDiagonal) 
 CYzTAYzCAYd = Matmul2(Conjg(Yz),TAYzCAYd,OnlyDiagonal) 
 CYzTAYzCAYs = Matmul2(Conjg(Yz),TAYzCAYs,OnlyDiagonal) 
 CYzTAYzCAYz = Matmul2(Conjg(Yz),TAYzCAYz,OnlyDiagonal) 
 CAYdTYdCYd = Matmul2(Conjg(AYd),TYdCYd,OnlyDiagonal) 
 CAYdTYdCYs = Matmul2(Conjg(AYd),TYdCYs,OnlyDiagonal) 
 CAYdTYdCYz = Matmul2(Conjg(AYd),TYdCYz,OnlyDiagonal) 
 CAYdTAYdCYd = Matmul2(Conjg(AYd),TAYdCYd,OnlyDiagonal) 
 CAYdTAYdCYs = Matmul2(Conjg(AYd),TAYdCYs,OnlyDiagonal) 
 CAYdTAYdCYz = Matmul2(Conjg(AYd),TAYdCYz,OnlyDiagonal) 
 CAYeTYeCYe = Matmul2(Conjg(AYe),TYeCYe,OnlyDiagonal) 
 CAYeTAYeCYe = Matmul2(Conjg(AYe),TAYeCYe,OnlyDiagonal) 
 CAYsYsCYd = Matmul2(adjAYs,YsCYd,OnlyDiagonal) 
 CAYsYsCYs = Matmul2(adjAYs,YsCYs,OnlyDiagonal) 
 CAYsYsCYz = Matmul2(adjAYs,YsCYz,OnlyDiagonal) 
 CAYsAYdadjYd = Matmul2(adjAYs,AYdadjYd,OnlyDiagonal) 
 CAYsAYsCYd = Matmul2(adjAYs,AYsCYd,OnlyDiagonal) 
 CAYsAYsCYs = Matmul2(adjAYs,AYsCYs,OnlyDiagonal) 
 CAYsAYsCYz = Matmul2(adjAYs,AYsCYz,OnlyDiagonal) 
 CAYsAYzadjYz = Matmul2(adjAYs,AYzadjYz,OnlyDiagonal) 
 CAYtYtadjYe = Matmul2(adjAYt,YtadjYe,OnlyDiagonal) 
 CAYtYtadjYz = Matmul2(adjAYt,YtadjYz,OnlyDiagonal) 
 CAYtYtCYt = Matmul2(adjAYt,YtCYt,OnlyDiagonal) 
 CAYtAYtadjYe = Matmul2(adjAYt,AYtadjYe,OnlyDiagonal) 
 CAYtAYtadjYz = Matmul2(adjAYt,AYtadjYz,OnlyDiagonal) 
 CAYtAYtCYt = Matmul2(adjAYt,AYtCYt,OnlyDiagonal) 
 CAYuTYuCYu = Matmul2(Conjg(AYu),TYuCYu,OnlyDiagonal) 
 CAYuTAYuCYu = Matmul2(Conjg(AYu),TAYuCYu,OnlyDiagonal) 
 CAYzTYzCYd = Matmul2(Conjg(AYz),TYzCYd,OnlyDiagonal) 
 CAYzTYzCYs = Matmul2(Conjg(AYz),TYzCYs,OnlyDiagonal) 
 CAYzTYzCYz = Matmul2(Conjg(AYz),TYzCYz,OnlyDiagonal) 
 CAYzTAYzCYd = Matmul2(Conjg(AYz),TAYzCYd,OnlyDiagonal) 
 CAYzTAYzCYs = Matmul2(Conjg(AYz),TAYzCYs,OnlyDiagonal) 
 CAYzTAYzCYz = Matmul2(Conjg(AYz),TAYzCYz,OnlyDiagonal) 
 TYdmd2CYs = Matmul2(Transpose(Yd),md2CYs,OnlyDiagonal) 
 TYdmd2CYz = Matmul2(Transpose(Yd),md2CYz,OnlyDiagonal) 
 TYdCYdTYd = Matmul2(Transpose(Yd),CYdTYd,OnlyDiagonal) 
 TYdCYdTAYd = Matmul2(Transpose(Yd),CYdTAYd,OnlyDiagonal) 
 TYdCYsmd2 = Matmul2(Transpose(Yd),CYsmd2,OnlyDiagonal) 
 TYdCYsYd = Matmul2(Transpose(Yd),CYsYd,OnlyDiagonal) 
 TYdCYsYs = Matmul2(Transpose(Yd),CYsYs,OnlyDiagonal) 
 TYdCYsYz = Matmul2(Transpose(Yd),CYsYz,OnlyDiagonal) 
 TYdCYsAYd = Matmul2(Transpose(Yd),CYsAYd,OnlyDiagonal) 
 TYdCYsAYs = Matmul2(Transpose(Yd),CYsAYs,OnlyDiagonal) 
 TYdCYsAYz = Matmul2(Transpose(Yd),CYsAYz,OnlyDiagonal) 
 TYdCYzml2 = Matmul2(Transpose(Yd),CYzml2,OnlyDiagonal) 
 TYdCYzYt = Matmul2(Transpose(Yd),CYzYt,OnlyDiagonal) 
 TYdCYzAYt = Matmul2(Transpose(Yd),CYzAYt,OnlyDiagonal) 
 TYeCYeTYz = Matmul2(Transpose(Ye),CYeTYz,OnlyDiagonal) 
 TYeCYeTAYz = Matmul2(Transpose(Ye),CYeTAYz,OnlyDiagonal) 
 TYuCYuTYd = Matmul2(Transpose(Yu),CYuTYd,OnlyDiagonal) 
 TYuCYuTAYd = Matmul2(Transpose(Yu),CYuTAYd,OnlyDiagonal) 
 TYzmd2CYd = Matmul2(Transpose(Yz),md2CYd,OnlyDiagonal) 
 TYzmd2CYs = Matmul2(Transpose(Yz),md2CYs,OnlyDiagonal) 
 TYzCYdmq2 = Matmul2(Transpose(Yz),CYdmq2,OnlyDiagonal) 
 TYzCYsmd2 = Matmul2(Transpose(Yz),CYsmd2,OnlyDiagonal) 
 TYzCYsYd = Matmul2(Transpose(Yz),CYsYd,OnlyDiagonal) 
 TYzCYsYs = Matmul2(Transpose(Yz),CYsYs,OnlyDiagonal) 
 TYzCYsYz = Matmul2(Transpose(Yz),CYsYz,OnlyDiagonal) 
 TYzCYsAYd = Matmul2(Transpose(Yz),CYsAYd,OnlyDiagonal) 
 TYzCYsAYs = Matmul2(Transpose(Yz),CYsAYs,OnlyDiagonal) 
 TYzCYsAYz = Matmul2(Transpose(Yz),CYsAYz,OnlyDiagonal) 
 TYzCYzTYz = Matmul2(Transpose(Yz),CYzTYz,OnlyDiagonal) 
 TYzCYzTAYz = Matmul2(Transpose(Yz),CYzTAYz,OnlyDiagonal) 
 TAYdCYdTYd = Matmul2(Transpose(AYd),CYdTYd,OnlyDiagonal) 
 TAYdCYsYd = Matmul2(Transpose(AYd),CYsYd,OnlyDiagonal) 
 TAYdCYsYs = Matmul2(Transpose(AYd),CYsYs,OnlyDiagonal) 
 TAYdCYsYz = Matmul2(Transpose(AYd),CYsYz,OnlyDiagonal) 
 TAYdCYzYt = Matmul2(Transpose(AYd),CYzYt,OnlyDiagonal) 
 TAYeCYeTYz = Matmul2(Transpose(AYe),CYeTYz,OnlyDiagonal) 
 TAYuCYuTYd = Matmul2(Transpose(AYu),CYuTYd,OnlyDiagonal) 
 TAYzCYsYd = Matmul2(Transpose(AYz),CYsYd,OnlyDiagonal) 
 TAYzCYsYs = Matmul2(Transpose(AYz),CYsYs,OnlyDiagonal) 
 TAYzCYsYz = Matmul2(Transpose(AYz),CYsYz,OnlyDiagonal) 
 TAYzCYzTYz = Matmul2(Transpose(AYz),CYzTYz,OnlyDiagonal) 
 md2YsCYsYs = Matmul2(md2,YsCYsYs,OnlyDiagonal) 
 md2CYdTYdCYd = Matmul2(md2,CYdTYdCYd,OnlyDiagonal) 
 md2CYdTYdCYs = Matmul2(md2,CYdTYdCYs,OnlyDiagonal) 
 md2CYdTYdCYz = Matmul2(md2,CYdTYdCYz,OnlyDiagonal) 
 md2CYsYdadjYd = Matmul2(md2,CYsYdadjYd,OnlyDiagonal) 
 md2CYsYsCYd = Matmul2(md2,CYsYsCYd,OnlyDiagonal) 
 md2CYsYsCYs = Matmul2(md2,CYsYsCYs,OnlyDiagonal) 
 md2CYsYsCYz = Matmul2(md2,CYsYsCYz,OnlyDiagonal) 
 md2CYsYzadjYz = Matmul2(md2,CYsYzadjYz,OnlyDiagonal) 
 md2CYzTYzCYd = Matmul2(md2,CYzTYzCYd,OnlyDiagonal) 
 md2CYzTYzCYs = Matmul2(md2,CYzTYzCYs,OnlyDiagonal) 
 md2CYzTYzCYz = Matmul2(md2,CYzTYzCYz,OnlyDiagonal) 
 me2CYeTYeCYe = Matmul2(me2,CYeTYeCYe,OnlyDiagonal) 
 ml2YtCYtYt = Matmul2(ml2,YtCYtYt,OnlyDiagonal) 
 ml2adjYeYeadjYe = Matmul2(ml2,adjYeYeadjYe,OnlyDiagonal) 
 ml2adjYeYeadjYz = Matmul2(ml2,adjYeYeadjYz,OnlyDiagonal) 
 ml2adjYeYeCYt = Matmul2(ml2,adjYeYeCYt,OnlyDiagonal) 
 ml2adjYzYdadjYd = Matmul2(ml2,adjYzYdadjYd,OnlyDiagonal) 
 ml2adjYzYsCYs = Matmul2(ml2,adjYzYsCYs,OnlyDiagonal) 
 ml2adjYzYzadjYe = Matmul2(ml2,adjYzYzadjYe,OnlyDiagonal) 
 ml2adjYzYzadjYz = Matmul2(ml2,adjYzYzadjYz,OnlyDiagonal) 
 ml2adjYzYzCYt = Matmul2(ml2,adjYzYzCYt,OnlyDiagonal) 
 ml2CYtYtadjYe = Matmul2(ml2,CYtYtadjYe,OnlyDiagonal) 
 ml2CYtYtadjYz = Matmul2(ml2,CYtYtadjYz,OnlyDiagonal) 
 ml2CYtYtCYt = Matmul2(ml2,CYtYtCYt,OnlyDiagonal) 
 mq2adjYdYdadjYd = Matmul2(mq2,adjYdYdadjYd,OnlyDiagonal) 
 mq2adjYdYdadjYu = Matmul2(mq2,adjYdYdadjYu,OnlyDiagonal) 
 mq2adjYdYsCYs = Matmul2(mq2,adjYdYsCYs,OnlyDiagonal) 
 mq2adjYdYzadjYz = Matmul2(mq2,adjYdYzadjYz,OnlyDiagonal) 
 mq2adjYuYuadjYd = Matmul2(mq2,adjYuYuadjYd,OnlyDiagonal) 
 mq2adjYuYuadjYu = Matmul2(mq2,adjYuYuadjYu,OnlyDiagonal) 
 mu2CYuTYuCYu = Matmul2(mu2,CYuTYuCYu,OnlyDiagonal) 
 Ydmq2adjYdYs = Matmul2(Yd,mq2adjYdYs,OnlyDiagonal) 
 YdadjYdmd2Ys = Matmul2(Yd,adjYdmd2Ys,OnlyDiagonal) 
 YdadjYdYdadjYd = Matmul2(Yd,adjYdYdadjYd,OnlyDiagonal) 
 YdadjYdYsmd2 = Matmul2(Yd,adjYdYsmd2,OnlyDiagonal) 
 YdadjYdYsCYs = Matmul2(Yd,adjYdYsCYs,OnlyDiagonal) 
 YdadjYdYzadjYz = Matmul2(Yd,adjYdYzadjYz,OnlyDiagonal) 
 YdadjYdAYdadjYd = Matmul2(Yd,adjYdAYdadjYd,OnlyDiagonal) 
 YdadjYdAYdadjAYd = Matmul2(Yd,adjYdAYdadjAYd,OnlyDiagonal) 
 YdadjYdAYsCYs = Matmul2(Yd,adjYdAYsCYs,OnlyDiagonal) 
 YdadjYdAYsCAYs = Matmul2(Yd,adjYdAYsCAYs,OnlyDiagonal) 
 YdadjYdAYzadjYz = Matmul2(Yd,adjYdAYzadjYz,OnlyDiagonal) 
 YdadjYdAYzadjAYz = Matmul2(Yd,adjYdAYzadjAYz,OnlyDiagonal) 
 YdadjYuYuadjYd = Matmul2(Yd,adjYuYuadjYd,OnlyDiagonal) 
 YdadjYuAYuadjYd = Matmul2(Yd,adjYuAYuadjYd,OnlyDiagonal) 
 YdadjYuAYuadjAYd = Matmul2(Yd,adjYuAYuadjAYd,OnlyDiagonal) 
 YdadjAYdAYdadjYd = Matmul2(Yd,adjAYdAYdadjYd,OnlyDiagonal) 
 YdadjAYdAYsCYs = Matmul2(Yd,adjAYdAYsCYs,OnlyDiagonal) 
 YdadjAYdAYzadjYz = Matmul2(Yd,adjAYdAYzadjYz,OnlyDiagonal) 
 YdadjAYuAYuadjYd = Matmul2(Yd,adjAYuAYuadjYd,OnlyDiagonal) 
 YeadjYeYeadjYe = Matmul2(Ye,adjYeYeadjYe,OnlyDiagonal) 
 YeadjYeAYeadjYe = Matmul2(Ye,adjYeAYeadjYe,OnlyDiagonal) 
 YeadjYeAYeadjAYe = Matmul2(Ye,adjYeAYeadjAYe,OnlyDiagonal) 
 YeadjYzYzadjYe = Matmul2(Ye,adjYzYzadjYe,OnlyDiagonal) 
 YeadjYzAYzadjYe = Matmul2(Ye,adjYzAYzadjYe,OnlyDiagonal) 
 YeadjYzAYzadjAYe = Matmul2(Ye,adjYzAYzadjAYe,OnlyDiagonal) 
 YeadjAYeAYeadjYe = Matmul2(Ye,adjAYeAYeadjYe,OnlyDiagonal) 
 YeadjAYzAYzadjYe = Matmul2(Ye,adjAYzAYzadjYe,OnlyDiagonal) 
 YeCYtYtadjYe = Matmul2(Ye,CYtYtadjYe,OnlyDiagonal) 
 YeCYtAYtadjYe = Matmul2(Ye,CYtAYtadjYe,OnlyDiagonal) 
 YeCYtAYtadjAYe = Matmul2(Ye,CYtAYtadjAYe,OnlyDiagonal) 
 YeCAYtAYtadjYe = Matmul2(Ye,CAYtAYtadjYe,OnlyDiagonal) 
 Ysmd2CYsYs = Matmul2(Ys,md2CYsYs,OnlyDiagonal) 
 YsCYdTYdCYs = Matmul2(Ys,CYdTYdCYs,OnlyDiagonal) 
 YsCYdTAYdCAYs = Matmul2(Ys,CYdTAYdCAYs,OnlyDiagonal) 
 YsCYsmd2Ys = Matmul2(Ys,CYsmd2Ys,OnlyDiagonal) 
 YsCYsYdadjYd = Matmul2(Ys,CYsYdadjYd,OnlyDiagonal) 
 YsCYsYsmd2 = Matmul2(Ys,CYsYsmd2,OnlyDiagonal) 
 YsCYsYsCYs = Matmul2(Ys,CYsYsCYs,OnlyDiagonal) 
 YsCYsYzadjYz = Matmul2(Ys,CYsYzadjYz,OnlyDiagonal) 
 YsCYsAYdadjYd = Matmul2(Ys,CYsAYdadjYd,OnlyDiagonal) 
 YsCYsAYdadjAYd = Matmul2(Ys,CYsAYdadjAYd,OnlyDiagonal) 
 YsCYsAYsCYs = Matmul2(Ys,CYsAYsCYs,OnlyDiagonal) 
 YsCYsAYsCAYs = Matmul2(Ys,CYsAYsCAYs,OnlyDiagonal) 
 YsCYsAYzadjYz = Matmul2(Ys,CYsAYzadjYz,OnlyDiagonal) 
 YsCYsAYzadjAYz = Matmul2(Ys,CYsAYzadjAYz,OnlyDiagonal) 
 YsCYzTYzCYs = Matmul2(Ys,CYzTYzCYs,OnlyDiagonal) 
 YsCYzTAYzCAYs = Matmul2(Ys,CYzTAYzCAYs,OnlyDiagonal) 
 YsCAYdTAYdCYs = Matmul2(Ys,CAYdTAYdCYs,OnlyDiagonal) 
 YsCAYsAYdadjYd = Matmul2(Ys,CAYsAYdadjYd,OnlyDiagonal) 
 YsCAYsAYsCYs = Matmul2(Ys,CAYsAYsCYs,OnlyDiagonal) 
 YsCAYsAYzadjYz = Matmul2(Ys,CAYsAYzadjYz,OnlyDiagonal) 
 YsCAYzTAYzCYs = Matmul2(Ys,CAYzTAYzCYs,OnlyDiagonal) 
 Ytml2adjYeYe = Matmul2(Yt,ml2adjYeYe,OnlyDiagonal) 
 Ytml2adjYzYz = Matmul2(Yt,ml2adjYzYz,OnlyDiagonal) 
 Ytml2CYtYt = Matmul2(Yt,ml2CYtYt,OnlyDiagonal) 
 YtadjYeme2Ye = Matmul2(Yt,adjYeme2Ye,OnlyDiagonal) 
 YtadjYeYeml2 = Matmul2(Yt,adjYeYeml2,OnlyDiagonal) 
 YtadjYeYeCYt = Matmul2(Yt,adjYeYeCYt,OnlyDiagonal) 
 YtadjYeAYeCYt = Matmul2(Yt,adjYeAYeCYt,OnlyDiagonal) 
 YtadjYeAYeCAYt = Matmul2(Yt,adjYeAYeCAYt,OnlyDiagonal) 
 YtadjYzmd2Yz = Matmul2(Yt,adjYzmd2Yz,OnlyDiagonal) 
 YtadjYzYzml2 = Matmul2(Yt,adjYzYzml2,OnlyDiagonal) 
 YtadjYzYzCYt = Matmul2(Yt,adjYzYzCYt,OnlyDiagonal) 
 YtadjYzAYzCYt = Matmul2(Yt,adjYzAYzCYt,OnlyDiagonal) 
 YtadjYzAYzCAYt = Matmul2(Yt,adjYzAYzCAYt,OnlyDiagonal) 
 YtadjAYeAYeCYt = Matmul2(Yt,adjAYeAYeCYt,OnlyDiagonal) 
 YtadjAYzAYzCYt = Matmul2(Yt,adjAYzAYzCYt,OnlyDiagonal) 
 YtCYtml2Yt = Matmul2(Yt,CYtml2Yt,OnlyDiagonal) 
 YtCYtYtml2 = Matmul2(Yt,CYtYtml2,OnlyDiagonal) 
 YtCYtYtCYt = Matmul2(Yt,CYtYtCYt,OnlyDiagonal) 
 YtCYtAYtCYt = Matmul2(Yt,CYtAYtCYt,OnlyDiagonal) 
 YtCYtAYtCAYt = Matmul2(Yt,CYtAYtCAYt,OnlyDiagonal) 
 YtCAYtAYtCYt = Matmul2(Yt,CAYtAYtCYt,OnlyDiagonal) 
 YuadjYdYdadjYu = Matmul2(Yu,adjYdYdadjYu,OnlyDiagonal) 
 YuadjYdAYdadjYu = Matmul2(Yu,adjYdAYdadjYu,OnlyDiagonal) 
 YuadjYdAYdadjAYu = Matmul2(Yu,adjYdAYdadjAYu,OnlyDiagonal) 
 YuadjYuYuadjYu = Matmul2(Yu,adjYuYuadjYu,OnlyDiagonal) 
 YuadjYuAYuadjYu = Matmul2(Yu,adjYuAYuadjYu,OnlyDiagonal) 
 YuadjYuAYuadjAYu = Matmul2(Yu,adjYuAYuadjAYu,OnlyDiagonal) 
 YuadjAYdAYdadjYu = Matmul2(Yu,adjAYdAYdadjYu,OnlyDiagonal) 
 YuadjAYuAYuadjYu = Matmul2(Yu,adjAYuAYuadjYu,OnlyDiagonal) 
 Yzml2adjYzYs = Matmul2(Yz,ml2adjYzYs,OnlyDiagonal) 
 YzadjYeYeadjYz = Matmul2(Yz,adjYeYeadjYz,OnlyDiagonal) 
 YzadjYeAYeadjYz = Matmul2(Yz,adjYeAYeadjYz,OnlyDiagonal) 
 YzadjYeAYeadjAYz = Matmul2(Yz,adjYeAYeadjAYz,OnlyDiagonal) 
 YzadjYzmd2Ys = Matmul2(Yz,adjYzmd2Ys,OnlyDiagonal) 
 YzadjYzYdadjYd = Matmul2(Yz,adjYzYdadjYd,OnlyDiagonal) 
 YzadjYzYsmd2 = Matmul2(Yz,adjYzYsmd2,OnlyDiagonal) 
 YzadjYzYsCYs = Matmul2(Yz,adjYzYsCYs,OnlyDiagonal) 
 YzadjYzYzadjYz = Matmul2(Yz,adjYzYzadjYz,OnlyDiagonal) 
 YzadjYzAYdadjYd = Matmul2(Yz,adjYzAYdadjYd,OnlyDiagonal) 
 YzadjYzAYdadjAYd = Matmul2(Yz,adjYzAYdadjAYd,OnlyDiagonal) 
 YzadjYzAYsCYs = Matmul2(Yz,adjYzAYsCYs,OnlyDiagonal) 
 YzadjYzAYsCAYs = Matmul2(Yz,adjYzAYsCAYs,OnlyDiagonal) 
 YzadjYzAYzadjYz = Matmul2(Yz,adjYzAYzadjYz,OnlyDiagonal) 
 YzadjYzAYzadjAYz = Matmul2(Yz,adjYzAYzadjAYz,OnlyDiagonal) 
 YzadjAYeAYeadjYz = Matmul2(Yz,adjAYeAYeadjYz,OnlyDiagonal) 
 YzadjAYzAYdadjYd = Matmul2(Yz,adjAYzAYdadjYd,OnlyDiagonal) 
 YzadjAYzAYsCYs = Matmul2(Yz,adjAYzAYsCYs,OnlyDiagonal) 
 YzadjAYzAYzadjYz = Matmul2(Yz,adjAYzAYzadjYz,OnlyDiagonal) 
 YzCYtYtadjYz = Matmul2(Yz,CYtYtadjYz,OnlyDiagonal) 
 YzCYtAYtadjYz = Matmul2(Yz,CYtAYtadjYz,OnlyDiagonal) 
 YzCYtAYtadjAYz = Matmul2(Yz,CYtAYtadjAYz,OnlyDiagonal) 
 YzCAYtAYtadjYz = Matmul2(Yz,CAYtAYtadjYz,OnlyDiagonal) 
 AYdadjYdYdadjAYd = Matmul2(AYd,adjYdYdadjAYd,OnlyDiagonal) 
 AYdadjYuYuadjAYd = Matmul2(AYd,adjYuYuadjAYd,OnlyDiagonal) 
 AYdadjAYdYdadjYd = Matmul2(AYd,adjAYdYdadjYd,OnlyDiagonal) 
 AYdadjAYuYuadjYd = Matmul2(AYd,adjAYuYuadjYd,OnlyDiagonal) 
 AYeadjYeYeadjAYe = Matmul2(AYe,adjYeYeadjAYe,OnlyDiagonal) 
 AYeadjYzYzadjAYe = Matmul2(AYe,adjYzYzadjAYe,OnlyDiagonal) 
 AYeadjAYeYeadjYe = Matmul2(AYe,adjAYeYeadjYe,OnlyDiagonal) 
 AYeadjAYzYzadjYe = Matmul2(AYe,adjAYzYzadjYe,OnlyDiagonal) 
 AYeCYtYtadjAYe = Matmul2(AYe,CYtYtadjAYe,OnlyDiagonal) 
 AYeCAYtYtadjYe = Matmul2(AYe,CAYtYtadjYe,OnlyDiagonal) 
 AYsCYdTYdCAYs = Matmul2(AYs,CYdTYdCAYs,OnlyDiagonal) 
 AYsCYsYsCAYs = Matmul2(AYs,CYsYsCAYs,OnlyDiagonal) 
 AYsCYzTYzCAYs = Matmul2(AYs,CYzTYzCAYs,OnlyDiagonal) 
 AYsCAYdTYdCYs = Matmul2(AYs,CAYdTYdCYs,OnlyDiagonal) 
 AYsCAYsYsCYs = Matmul2(AYs,CAYsYsCYs,OnlyDiagonal) 
 AYsCAYzTYzCYs = Matmul2(AYs,CAYzTYzCYs,OnlyDiagonal) 
 AYtadjYeYeCAYt = Matmul2(AYt,adjYeYeCAYt,OnlyDiagonal) 
 AYtadjYzYzCAYt = Matmul2(AYt,adjYzYzCAYt,OnlyDiagonal) 
 AYtadjAYeYeCYt = Matmul2(AYt,adjAYeYeCYt,OnlyDiagonal) 
 AYtadjAYzYzCYt = Matmul2(AYt,adjAYzYzCYt,OnlyDiagonal) 
 AYtCYtYtCAYt = Matmul2(AYt,CYtYtCAYt,OnlyDiagonal) 
 AYtCAYtYtCYt = Matmul2(AYt,CAYtYtCYt,OnlyDiagonal) 
 AYuadjYdYdadjAYu = Matmul2(AYu,adjYdYdadjAYu,OnlyDiagonal) 
 AYuadjYuYuadjAYu = Matmul2(AYu,adjYuYuadjAYu,OnlyDiagonal) 
 AYuadjAYdYdadjYu = Matmul2(AYu,adjAYdYdadjYu,OnlyDiagonal) 
 AYuadjAYuYuadjYu = Matmul2(AYu,adjAYuYuadjYu,OnlyDiagonal) 
 AYzadjYeYeadjAYz = Matmul2(AYz,adjYeYeadjAYz,OnlyDiagonal) 
 AYzadjYzYzadjAYz = Matmul2(AYz,adjYzYzadjAYz,OnlyDiagonal) 
 AYzadjAYeYeadjYz = Matmul2(AYz,adjAYeYeadjYz,OnlyDiagonal) 
 AYzadjAYzYzadjYz = Matmul2(AYz,adjAYzYzadjYz,OnlyDiagonal) 
 AYzCYtYtadjAYz = Matmul2(AYz,CYtYtadjAYz,OnlyDiagonal) 
 AYzCAYtYtadjYz = Matmul2(AYz,CAYtYtadjYz,OnlyDiagonal) 
 adjYdmd2YdadjYd = Matmul2(adjYd,md2YdadjYd,OnlyDiagonal) 
 adjYdmd2YdadjYu = Matmul2(adjYd,md2YdadjYu,OnlyDiagonal) 
 adjYdmd2YsCYs = Matmul2(adjYd,md2YsCYs,OnlyDiagonal) 
 adjYdmd2YzadjYz = Matmul2(adjYd,md2YzadjYz,OnlyDiagonal) 
 adjYdYdmq2adjYd = Matmul2(adjYd,Ydmq2adjYd,OnlyDiagonal) 
 adjYdYdmq2adjYu = Matmul2(adjYd,Ydmq2adjYu,OnlyDiagonal) 
 adjYdYdadjYdmd2 = Matmul2(adjYd,YdadjYdmd2,OnlyDiagonal) 
 adjYdYdadjYdYd = Matmul2(adjYd,YdadjYdYd,OnlyDiagonal) 
 adjYdYdadjYdYs = Matmul2(adjYd,YdadjYdYs,OnlyDiagonal) 
 adjYdYdadjYdYz = Matmul2(adjYd,YdadjYdYz,OnlyDiagonal) 
 adjYdYdadjYdAYd = Matmul2(adjYd,YdadjYdAYd,OnlyDiagonal) 
 adjYdYdadjYdAYs = Matmul2(adjYd,YdadjYdAYs,OnlyDiagonal) 
 adjYdYdadjYdAYz = Matmul2(adjYd,YdadjYdAYz,OnlyDiagonal) 
 adjYdYdadjYumu2 = Matmul2(adjYd,YdadjYumu2,OnlyDiagonal) 
 adjYdYdadjYuYu = Matmul2(adjYd,YdadjYuYu,OnlyDiagonal) 
 adjYdYdadjYuAYu = Matmul2(adjYd,YdadjYuAYu,OnlyDiagonal) 
 adjYdYsCYsYd = Matmul2(adjYd,YsCYsYd,OnlyDiagonal) 
 adjYdYsCYsAYd = Matmul2(adjYd,YsCYsAYd,OnlyDiagonal) 
 adjYdYzml2adjYz = Matmul2(adjYd,Yzml2adjYz,OnlyDiagonal) 
 adjYdYzadjYzYd = Matmul2(adjYd,YzadjYzYd,OnlyDiagonal) 
 adjYdYzadjYzAYd = Matmul2(adjYd,YzadjYzAYd,OnlyDiagonal) 
 adjYdAYdadjYdYd = Matmul2(adjYd,AYdadjYdYd,OnlyDiagonal) 
 adjYdAYdadjYdYs = Matmul2(adjYd,AYdadjYdYs,OnlyDiagonal) 
 adjYdAYdadjYdYz = Matmul2(adjYd,AYdadjYdYz,OnlyDiagonal) 
 adjYdAYdadjYuYu = Matmul2(adjYd,AYdadjYuYu,OnlyDiagonal) 
 adjYdAYsCYsYd = Matmul2(adjYd,AYsCYsYd,OnlyDiagonal) 
 adjYdAYzadjYzYd = Matmul2(adjYd,AYzadjYzYd,OnlyDiagonal) 
 adjYeme2YeadjYe = Matmul2(adjYe,me2YeadjYe,OnlyDiagonal) 
 adjYeme2YeadjYz = Matmul2(adjYe,me2YeadjYz,OnlyDiagonal) 
 adjYeme2YeCYt = Matmul2(adjYe,me2YeCYt,OnlyDiagonal) 
 adjYeYeml2adjYe = Matmul2(adjYe,Yeml2adjYe,OnlyDiagonal) 
 adjYeYeml2adjYz = Matmul2(adjYe,Yeml2adjYz,OnlyDiagonal) 
 adjYeYeml2CYt = Matmul2(adjYe,Yeml2CYt,OnlyDiagonal) 
 adjYeYeadjYeme2 = Matmul2(adjYe,YeadjYeme2,OnlyDiagonal) 
 adjYeYeadjYeYe = Matmul2(adjYe,YeadjYeYe,OnlyDiagonal) 
 adjYeYeadjYeAYe = Matmul2(adjYe,YeadjYeAYe,OnlyDiagonal) 
 adjYeYeadjYzmd2 = Matmul2(adjYe,YeadjYzmd2,OnlyDiagonal) 
 adjYeYeadjYzYd = Matmul2(adjYe,YeadjYzYd,OnlyDiagonal) 
 adjYeYeadjYzYs = Matmul2(adjYe,YeadjYzYs,OnlyDiagonal) 
 adjYeYeadjYzYz = Matmul2(adjYe,YeadjYzYz,OnlyDiagonal) 
 adjYeYeadjYzAYd = Matmul2(adjYe,YeadjYzAYd,OnlyDiagonal) 
 adjYeYeadjYzAYs = Matmul2(adjYe,YeadjYzAYs,OnlyDiagonal) 
 adjYeYeadjYzAYz = Matmul2(adjYe,YeadjYzAYz,OnlyDiagonal) 
 adjYeYeCYtml2 = Matmul2(adjYe,YeCYtml2,OnlyDiagonal) 
 adjYeYeCYtYt = Matmul2(adjYe,YeCYtYt,OnlyDiagonal) 
 adjYeYeCYtAYt = Matmul2(adjYe,YeCYtAYt,OnlyDiagonal) 
 adjYeAYeadjYeYe = Matmul2(adjYe,AYeadjYeYe,OnlyDiagonal) 
 adjYeAYeadjYzYd = Matmul2(adjYe,AYeadjYzYd,OnlyDiagonal) 
 adjYeAYeadjYzYs = Matmul2(adjYe,AYeadjYzYs,OnlyDiagonal) 
 adjYeAYeadjYzYz = Matmul2(adjYe,AYeadjYzYz,OnlyDiagonal) 
 adjYeAYeCYtYt = Matmul2(adjYe,AYeCYtYt,OnlyDiagonal) 
 adjYumu2YuadjYd = Matmul2(adjYu,mu2YuadjYd,OnlyDiagonal) 
 adjYumu2YuadjYu = Matmul2(adjYu,mu2YuadjYu,OnlyDiagonal) 
 adjYuYumq2adjYd = Matmul2(adjYu,Yumq2adjYd,OnlyDiagonal) 
 adjYuYumq2adjYu = Matmul2(adjYu,Yumq2adjYu,OnlyDiagonal) 
 adjYuYuadjYdmd2 = Matmul2(adjYu,YuadjYdmd2,OnlyDiagonal) 
 adjYuYuadjYdYd = Matmul2(adjYu,YuadjYdYd,OnlyDiagonal) 
 adjYuYuadjYdYs = Matmul2(adjYu,YuadjYdYs,OnlyDiagonal) 
 adjYuYuadjYdYz = Matmul2(adjYu,YuadjYdYz,OnlyDiagonal) 
 adjYuYuadjYdAYd = Matmul2(adjYu,YuadjYdAYd,OnlyDiagonal) 
 adjYuYuadjYdAYs = Matmul2(adjYu,YuadjYdAYs,OnlyDiagonal) 
 adjYuYuadjYdAYz = Matmul2(adjYu,YuadjYdAYz,OnlyDiagonal) 
 adjYuYuadjYumu2 = Matmul2(adjYu,YuadjYumu2,OnlyDiagonal) 
 adjYuYuadjYuYu = Matmul2(adjYu,YuadjYuYu,OnlyDiagonal) 
 adjYuYuadjYuAYu = Matmul2(adjYu,YuadjYuAYu,OnlyDiagonal) 
 adjYuAYuadjYdYd = Matmul2(adjYu,AYuadjYdYd,OnlyDiagonal) 
 adjYuAYuadjYdYs = Matmul2(adjYu,AYuadjYdYs,OnlyDiagonal) 
 adjYuAYuadjYdYz = Matmul2(adjYu,AYuadjYdYz,OnlyDiagonal) 
 adjYuAYuadjYuYu = Matmul2(adjYu,AYuadjYuYu,OnlyDiagonal) 
 adjYzmd2YdadjYd = Matmul2(adjYz,md2YdadjYd,OnlyDiagonal) 
 adjYzmd2YsCYs = Matmul2(adjYz,md2YsCYs,OnlyDiagonal) 
 adjYzmd2YzadjYe = Matmul2(adjYz,md2YzadjYe,OnlyDiagonal) 
 adjYzmd2YzadjYz = Matmul2(adjYz,md2YzadjYz,OnlyDiagonal) 
 adjYzmd2YzCYt = Matmul2(adjYz,md2YzCYt,OnlyDiagonal) 
 adjYzYdmq2adjYd = Matmul2(adjYz,Ydmq2adjYd,OnlyDiagonal) 
 adjYzYdadjYdYz = Matmul2(adjYz,YdadjYdYz,OnlyDiagonal) 
 adjYzYdadjYdAYz = Matmul2(adjYz,YdadjYdAYz,OnlyDiagonal) 
 adjYzYsCYsYz = Matmul2(adjYz,YsCYsYz,OnlyDiagonal) 
 adjYzYsCYsAYz = Matmul2(adjYz,YsCYsAYz,OnlyDiagonal) 
 adjYzYzml2adjYe = Matmul2(adjYz,Yzml2adjYe,OnlyDiagonal) 
 adjYzYzml2adjYz = Matmul2(adjYz,Yzml2adjYz,OnlyDiagonal) 
 adjYzYzml2CYt = Matmul2(adjYz,Yzml2CYt,OnlyDiagonal) 
 adjYzYzadjYeme2 = Matmul2(adjYz,YzadjYeme2,OnlyDiagonal) 
 adjYzYzadjYeYe = Matmul2(adjYz,YzadjYeYe,OnlyDiagonal) 
 adjYzYzadjYeAYe = Matmul2(adjYz,YzadjYeAYe,OnlyDiagonal) 
 adjYzYzadjYzmd2 = Matmul2(adjYz,YzadjYzmd2,OnlyDiagonal) 
 adjYzYzadjYzYd = Matmul2(adjYz,YzadjYzYd,OnlyDiagonal) 
 adjYzYzadjYzYs = Matmul2(adjYz,YzadjYzYs,OnlyDiagonal) 
 adjYzYzadjYzYz = Matmul2(adjYz,YzadjYzYz,OnlyDiagonal) 
 adjYzYzadjYzAYd = Matmul2(adjYz,YzadjYzAYd,OnlyDiagonal) 
 adjYzYzadjYzAYs = Matmul2(adjYz,YzadjYzAYs,OnlyDiagonal) 
 adjYzYzadjYzAYz = Matmul2(adjYz,YzadjYzAYz,OnlyDiagonal) 
 adjYzYzCYtml2 = Matmul2(adjYz,YzCYtml2,OnlyDiagonal) 
 adjYzYzCYtYt = Matmul2(adjYz,YzCYtYt,OnlyDiagonal) 
 adjYzYzCYtAYt = Matmul2(adjYz,YzCYtAYt,OnlyDiagonal) 
 adjYzAYdadjYdYz = Matmul2(adjYz,AYdadjYdYz,OnlyDiagonal) 
 adjYzAYsCYsYz = Matmul2(adjYz,AYsCYsYz,OnlyDiagonal) 
 adjYzAYzadjYeYe = Matmul2(adjYz,AYzadjYeYe,OnlyDiagonal) 
 adjYzAYzadjYzYd = Matmul2(adjYz,AYzadjYzYd,OnlyDiagonal) 
 adjYzAYzadjYzYs = Matmul2(adjYz,AYzadjYzYs,OnlyDiagonal) 
 adjYzAYzadjYzYz = Matmul2(adjYz,AYzadjYzYz,OnlyDiagonal) 
 adjYzAYzCYtYt = Matmul2(adjYz,AYzCYtYt,OnlyDiagonal) 
 CYdmq2TYdCYd = Matmul2(Conjg(Yd),mq2TYdCYd,OnlyDiagonal) 
 CYdmq2TYdCYs = Matmul2(Conjg(Yd),mq2TYdCYs,OnlyDiagonal) 
 CYdmq2TYdCYz = Matmul2(Conjg(Yd),mq2TYdCYz,OnlyDiagonal) 
 CYdTYdmd2CYd = Matmul2(Conjg(Yd),TYdmd2CYd,OnlyDiagonal) 
 CYdTYdmd2CYs = Matmul2(Conjg(Yd),TYdmd2CYs,OnlyDiagonal) 
 CYdTYdmd2CYz = Matmul2(Conjg(Yd),TYdmd2CYz,OnlyDiagonal) 
 CYdTYdCYdmq2 = Matmul2(Conjg(Yd),TYdCYdmq2,OnlyDiagonal) 
 CYdTYdCYdTYd = Matmul2(Conjg(Yd),TYdCYdTYd,OnlyDiagonal) 
 CYdTYdCYdTAYd = Matmul2(Conjg(Yd),TYdCYdTAYd,OnlyDiagonal) 
 CYdTYdCYsmd2 = Matmul2(Conjg(Yd),TYdCYsmd2,OnlyDiagonal) 
 CYdTYdCYsYd = Matmul2(Conjg(Yd),TYdCYsYd,OnlyDiagonal) 
 CYdTYdCYsYs = Matmul2(Conjg(Yd),TYdCYsYs,OnlyDiagonal) 
 CYdTYdCYsYz = Matmul2(Conjg(Yd),TYdCYsYz,OnlyDiagonal) 
 CYdTYdCYsAYd = Matmul2(Conjg(Yd),TYdCYsAYd,OnlyDiagonal) 
 CYdTYdCYsAYs = Matmul2(Conjg(Yd),TYdCYsAYs,OnlyDiagonal) 
 CYdTYdCYsAYz = Matmul2(Conjg(Yd),TYdCYsAYz,OnlyDiagonal) 
 CYdTYdCYzml2 = Matmul2(Conjg(Yd),TYdCYzml2,OnlyDiagonal) 
 CYdTYdCYzYt = Matmul2(Conjg(Yd),TYdCYzYt,OnlyDiagonal) 
 CYdTYdCYzAYt = Matmul2(Conjg(Yd),TYdCYzAYt,OnlyDiagonal) 
 CYdTYuCYuTYd = Matmul2(Conjg(Yd),TYuCYuTYd,OnlyDiagonal) 
 CYdTYuCYuTAYd = Matmul2(Conjg(Yd),TYuCYuTAYd,OnlyDiagonal) 
 CYdTAYdCYdTYd = Matmul2(Conjg(Yd),TAYdCYdTYd,OnlyDiagonal) 
 CYdTAYdCYsYd = Matmul2(Conjg(Yd),TAYdCYsYd,OnlyDiagonal) 
 CYdTAYdCYsYs = Matmul2(Conjg(Yd),TAYdCYsYs,OnlyDiagonal) 
 CYdTAYdCYsYz = Matmul2(Conjg(Yd),TAYdCYsYz,OnlyDiagonal) 
 CYdTAYdCYzYt = Matmul2(Conjg(Yd),TAYdCYzYt,OnlyDiagonal) 
 CYdTAYuCYuTYd = Matmul2(Conjg(Yd),TAYuCYuTYd,OnlyDiagonal) 
 CYeml2TYeCYe = Matmul2(Conjg(Ye),ml2TYeCYe,OnlyDiagonal) 
 CYeTYeme2CYe = Matmul2(Conjg(Ye),TYeme2CYe,OnlyDiagonal) 
 CYeTYeCYeml2 = Matmul2(Conjg(Ye),TYeCYeml2,OnlyDiagonal) 
 CYeTYeCYeYt = Matmul2(Conjg(Ye),TYeCYeYt,OnlyDiagonal) 
 CYeTYeCYeAYt = Matmul2(Conjg(Ye),TYeCYeAYt,OnlyDiagonal) 
 CYeTAYeCYeYt = Matmul2(Conjg(Ye),TAYeCYeYt,OnlyDiagonal) 
 CYsmd2YdadjYd = Matmul2(adjYs,md2YdadjYd,OnlyDiagonal) 
 CYsmd2YsCYd = Matmul2(adjYs,md2YsCYd,OnlyDiagonal) 
 CYsmd2YsCYs = Matmul2(adjYs,md2YsCYs,OnlyDiagonal) 
 CYsmd2YsCYz = Matmul2(adjYs,md2YsCYz,OnlyDiagonal) 
 CYsmd2YzadjYz = Matmul2(adjYs,md2YzadjYz,OnlyDiagonal) 
 CYsYdmq2adjYd = Matmul2(adjYs,Ydmq2adjYd,OnlyDiagonal) 
 CYsYdadjYdmd2 = Matmul2(adjYs,YdadjYdmd2,OnlyDiagonal) 
 CYsYdadjYdYs = Matmul2(adjYs,YdadjYdYs,OnlyDiagonal) 
 CYsYdadjYdAYs = Matmul2(adjYs,YdadjYdAYs,OnlyDiagonal) 
 CYsYdadjAYdAYs = Matmul2(adjYs,YdadjAYdAYs,OnlyDiagonal) 
 CYsYsmd2CYd = Matmul2(adjYs,Ysmd2CYd,OnlyDiagonal) 
 CYsYsmd2CYs = Matmul2(adjYs,Ysmd2CYs,OnlyDiagonal) 
 CYsYsmd2CYz = Matmul2(adjYs,Ysmd2CYz,OnlyDiagonal) 
 CYsYsCYdmq2 = Matmul2(adjYs,YsCYdmq2,OnlyDiagonal) 
 CYsYsCYsmd2 = Matmul2(adjYs,YsCYsmd2,OnlyDiagonal) 
 CYsYsCYsYd = Matmul2(adjYs,YsCYsYd,OnlyDiagonal) 
 CYsYsCYsYs = Matmul2(adjYs,YsCYsYs,OnlyDiagonal) 
 CYsYsCYsYz = Matmul2(adjYs,YsCYsYz,OnlyDiagonal) 
 CYsYsCYsAYd = Matmul2(adjYs,YsCYsAYd,OnlyDiagonal) 
 CYsYsCYsAYs = Matmul2(adjYs,YsCYsAYs,OnlyDiagonal) 
 CYsYsCYsAYz = Matmul2(adjYs,YsCYsAYz,OnlyDiagonal) 
 CYsYsCYzml2 = Matmul2(adjYs,YsCYzml2,OnlyDiagonal) 
 CYsYsCYzYt = Matmul2(adjYs,YsCYzYt,OnlyDiagonal) 
 CYsYsCYzAYt = Matmul2(adjYs,YsCYzAYt,OnlyDiagonal) 
 CYsYsCAYsAYs = Matmul2(adjYs,YsCAYsAYs,OnlyDiagonal) 
 CYsYzml2adjYz = Matmul2(adjYs,Yzml2adjYz,OnlyDiagonal) 
 CYsYzadjYzmd2 = Matmul2(adjYs,YzadjYzmd2,OnlyDiagonal) 
 CYsYzadjYzYs = Matmul2(adjYs,YzadjYzYs,OnlyDiagonal) 
 CYsYzadjYzAYs = Matmul2(adjYs,YzadjYzAYs,OnlyDiagonal) 
 CYsYzadjAYzAYs = Matmul2(adjYs,YzadjAYzAYs,OnlyDiagonal) 
 CYsAYdadjYdYs = Matmul2(adjYs,AYdadjYdYs,OnlyDiagonal) 
 CYsAYdadjAYdYs = Matmul2(adjYs,AYdadjAYdYs,OnlyDiagonal) 
 CYsAYsCYsYd = Matmul2(adjYs,AYsCYsYd,OnlyDiagonal) 
 CYsAYsCYsYs = Matmul2(adjYs,AYsCYsYs,OnlyDiagonal) 
 CYsAYsCYsYz = Matmul2(adjYs,AYsCYsYz,OnlyDiagonal) 
 CYsAYsCYzYt = Matmul2(adjYs,AYsCYzYt,OnlyDiagonal) 
 CYsAYsCAYsYs = Matmul2(adjYs,AYsCAYsYs,OnlyDiagonal) 
 CYsAYzadjYzYs = Matmul2(adjYs,AYzadjYzYs,OnlyDiagonal) 
 CYsAYzadjAYzYs = Matmul2(adjYs,AYzadjAYzYs,OnlyDiagonal) 
 CYtml2YtadjYe = Matmul2(adjYt,ml2YtadjYe,OnlyDiagonal) 
 CYtml2YtadjYz = Matmul2(adjYt,ml2YtadjYz,OnlyDiagonal) 
 CYtml2YtCYt = Matmul2(adjYt,ml2YtCYt,OnlyDiagonal) 
 CYtYtml2adjYe = Matmul2(adjYt,Ytml2adjYe,OnlyDiagonal) 
 CYtYtml2adjYz = Matmul2(adjYt,Ytml2adjYz,OnlyDiagonal) 
 CYtYtml2CYt = Matmul2(adjYt,Ytml2CYt,OnlyDiagonal) 
 CYtYtadjYeme2 = Matmul2(adjYt,YtadjYeme2,OnlyDiagonal) 
 CYtYtadjYeYe = Matmul2(adjYt,YtadjYeYe,OnlyDiagonal) 
 CYtYtadjYeAYe = Matmul2(adjYt,YtadjYeAYe,OnlyDiagonal) 
 CYtYtadjYzmd2 = Matmul2(adjYt,YtadjYzmd2,OnlyDiagonal) 
 CYtYtadjYzYd = Matmul2(adjYt,YtadjYzYd,OnlyDiagonal) 
 CYtYtadjYzYs = Matmul2(adjYt,YtadjYzYs,OnlyDiagonal) 
 CYtYtadjYzYz = Matmul2(adjYt,YtadjYzYz,OnlyDiagonal) 
 CYtYtadjYzAYd = Matmul2(adjYt,YtadjYzAYd,OnlyDiagonal) 
 CYtYtadjYzAYs = Matmul2(adjYt,YtadjYzAYs,OnlyDiagonal) 
 CYtYtadjYzAYz = Matmul2(adjYt,YtadjYzAYz,OnlyDiagonal) 
 CYtYtadjAYeAYe = Matmul2(adjYt,YtadjAYeAYe,OnlyDiagonal) 
 CYtYtadjAYzAYz = Matmul2(adjYt,YtadjAYzAYz,OnlyDiagonal) 
 CYtYtCYtml2 = Matmul2(adjYt,YtCYtml2,OnlyDiagonal) 
 CYtYtCYtYt = Matmul2(adjYt,YtCYtYt,OnlyDiagonal) 
 CYtYtCYtAYt = Matmul2(adjYt,YtCYtAYt,OnlyDiagonal) 
 CYtYtCAYtAYt = Matmul2(adjYt,YtCAYtAYt,OnlyDiagonal) 
 CYtAYtadjYeYe = Matmul2(adjYt,AYtadjYeYe,OnlyDiagonal) 
 CYtAYtadjYzYd = Matmul2(adjYt,AYtadjYzYd,OnlyDiagonal) 
 CYtAYtadjYzYs = Matmul2(adjYt,AYtadjYzYs,OnlyDiagonal) 
 CYtAYtadjYzYz = Matmul2(adjYt,AYtadjYzYz,OnlyDiagonal) 
 CYtAYtadjAYeYe = Matmul2(adjYt,AYtadjAYeYe,OnlyDiagonal) 
 CYtAYtadjAYzYz = Matmul2(adjYt,AYtadjAYzYz,OnlyDiagonal) 
 CYtAYtCYtYt = Matmul2(adjYt,AYtCYtYt,OnlyDiagonal) 
 CYtAYtCAYtYt = Matmul2(adjYt,AYtCAYtYt,OnlyDiagonal) 
 CYtTYeCYeYt = Matmul2(adjYt,TYeCYeYt,OnlyDiagonal) 
 CYtTYeCYeAYt = Matmul2(adjYt,TYeCYeAYt,OnlyDiagonal) 
 CYtTYzCYzYt = Matmul2(adjYt,TYzCYzYt,OnlyDiagonal) 
 CYtTYzCYzAYt = Matmul2(adjYt,TYzCYzAYt,OnlyDiagonal) 
 CYtTAYeCYeYt = Matmul2(adjYt,TAYeCYeYt,OnlyDiagonal) 
 CYtTAYzCYzYt = Matmul2(adjYt,TAYzCYzYt,OnlyDiagonal) 
 CYumq2TYuCYu = Matmul2(Conjg(Yu),mq2TYuCYu,OnlyDiagonal) 
 CYuTYumu2CYu = Matmul2(Conjg(Yu),TYumu2CYu,OnlyDiagonal) 
 CYuTYuCYumq2 = Matmul2(Conjg(Yu),TYuCYumq2,OnlyDiagonal) 
 CYzml2TYzCYd = Matmul2(Conjg(Yz),ml2TYzCYd,OnlyDiagonal) 
 CYzml2TYzCYs = Matmul2(Conjg(Yz),ml2TYzCYs,OnlyDiagonal) 
 CYzml2TYzCYz = Matmul2(Conjg(Yz),ml2TYzCYz,OnlyDiagonal) 
 CYzYtCYtTYz = Matmul2(Conjg(Yz),YtCYtTYz,OnlyDiagonal) 
 CYzYtCYtTAYz = Matmul2(Conjg(Yz),YtCYtTAYz,OnlyDiagonal) 
 CYzAYtCYtTYz = Matmul2(Conjg(Yz),AYtCYtTYz,OnlyDiagonal) 
 CYzTYeCYeTYz = Matmul2(Conjg(Yz),TYeCYeTYz,OnlyDiagonal) 
 CYzTYeCYeTAYz = Matmul2(Conjg(Yz),TYeCYeTAYz,OnlyDiagonal) 
 CYzTYzmd2CYd = Matmul2(Conjg(Yz),TYzmd2CYd,OnlyDiagonal) 
 CYzTYzmd2CYs = Matmul2(Conjg(Yz),TYzmd2CYs,OnlyDiagonal) 
 CYzTYzmd2CYz = Matmul2(Conjg(Yz),TYzmd2CYz,OnlyDiagonal) 
 CYzTYzCYdmq2 = Matmul2(Conjg(Yz),TYzCYdmq2,OnlyDiagonal) 
 CYzTYzCYsmd2 = Matmul2(Conjg(Yz),TYzCYsmd2,OnlyDiagonal) 
 CYzTYzCYsYd = Matmul2(Conjg(Yz),TYzCYsYd,OnlyDiagonal) 
 CYzTYzCYsYs = Matmul2(Conjg(Yz),TYzCYsYs,OnlyDiagonal) 
 CYzTYzCYsYz = Matmul2(Conjg(Yz),TYzCYsYz,OnlyDiagonal) 
 CYzTYzCYsAYd = Matmul2(Conjg(Yz),TYzCYsAYd,OnlyDiagonal) 
 CYzTYzCYsAYs = Matmul2(Conjg(Yz),TYzCYsAYs,OnlyDiagonal) 
 CYzTYzCYsAYz = Matmul2(Conjg(Yz),TYzCYsAYz,OnlyDiagonal) 
 CYzTYzCYzml2 = Matmul2(Conjg(Yz),TYzCYzml2,OnlyDiagonal) 
 CYzTYzCYzYt = Matmul2(Conjg(Yz),TYzCYzYt,OnlyDiagonal) 
 CYzTYzCYzAYt = Matmul2(Conjg(Yz),TYzCYzAYt,OnlyDiagonal) 
 CYzTYzCYzTYz = Matmul2(Conjg(Yz),TYzCYzTYz,OnlyDiagonal) 
 CYzTYzCYzTAYz = Matmul2(Conjg(Yz),TYzCYzTAYz,OnlyDiagonal) 
 CYzTAYeCYeTYz = Matmul2(Conjg(Yz),TAYeCYeTYz,OnlyDiagonal) 
 CYzTAYzCYsYd = Matmul2(Conjg(Yz),TAYzCYsYd,OnlyDiagonal) 
 CYzTAYzCYsYs = Matmul2(Conjg(Yz),TAYzCYsYs,OnlyDiagonal) 
 CYzTAYzCYsYz = Matmul2(Conjg(Yz),TAYzCYsYz,OnlyDiagonal) 
 CYzTAYzCYzYt = Matmul2(Conjg(Yz),TAYzCYzYt,OnlyDiagonal) 
 CYzTAYzCYzTYz = Matmul2(Conjg(Yz),TAYzCYzTYz,OnlyDiagonal) 
 TYdCYdTYdCYd = Matmul2(Transpose(Yd),CYdTYdCYd,OnlyDiagonal) 
 TYdCYdTAYdCAYd = Matmul2(Transpose(Yd),CYdTAYdCAYd,OnlyDiagonal) 
 TYdCYsYsCYd = Matmul2(Transpose(Yd),CYsYsCYd,OnlyDiagonal) 
 TYdCYsAYsCAYd = Matmul2(Transpose(Yd),CYsAYsCAYd,OnlyDiagonal) 
 TYdCYzTYzCYd = Matmul2(Transpose(Yd),CYzTYzCYd,OnlyDiagonal) 
 TYdCYzTAYzCAYd = Matmul2(Transpose(Yd),CYzTAYzCAYd,OnlyDiagonal) 
 TYdCAYdTAYdCYd = Matmul2(Transpose(Yd),CAYdTAYdCYd,OnlyDiagonal) 
 TYdCAYsAYsCYd = Matmul2(Transpose(Yd),CAYsAYsCYd,OnlyDiagonal) 
 TYdCAYzTAYzCYd = Matmul2(Transpose(Yd),CAYzTAYzCYd,OnlyDiagonal) 
 TYeCYeTYeCYe = Matmul2(Transpose(Ye),CYeTYeCYe,OnlyDiagonal) 
 TYeCYeTAYeCAYe = Matmul2(Transpose(Ye),CYeTAYeCAYe,OnlyDiagonal) 
 TYeCAYeTAYeCYe = Matmul2(Transpose(Ye),CAYeTAYeCYe,OnlyDiagonal) 
 TYuCYuTYuCYu = Matmul2(Transpose(Yu),CYuTYuCYu,OnlyDiagonal) 
 TYuCYuTAYuCAYu = Matmul2(Transpose(Yu),CYuTAYuCAYu,OnlyDiagonal) 
 TYuCAYuTAYuCYu = Matmul2(Transpose(Yu),CAYuTAYuCYu,OnlyDiagonal) 
 TYzCYdTYdCYz = Matmul2(Transpose(Yz),CYdTYdCYz,OnlyDiagonal) 
 TYzCYdTAYdCAYz = Matmul2(Transpose(Yz),CYdTAYdCAYz,OnlyDiagonal) 
 TYzCYsYsCYz = Matmul2(Transpose(Yz),CYsYsCYz,OnlyDiagonal) 
 TYzCYsAYsCAYz = Matmul2(Transpose(Yz),CYsAYsCAYz,OnlyDiagonal) 
 TYzCYzTYzCYz = Matmul2(Transpose(Yz),CYzTYzCYz,OnlyDiagonal) 
 TYzCYzTAYzCAYz = Matmul2(Transpose(Yz),CYzTAYzCAYz,OnlyDiagonal) 
 TYzCAYdTAYdCYz = Matmul2(Transpose(Yz),CAYdTAYdCYz,OnlyDiagonal) 
 TYzCAYsAYsCYz = Matmul2(Transpose(Yz),CAYsAYsCYz,OnlyDiagonal) 
 TYzCAYzTAYzCYz = Matmul2(Transpose(Yz),CAYzTAYzCYz,OnlyDiagonal) 
 TAYdCYdTYdCAYd = Matmul2(Transpose(AYd),CYdTYdCAYd,OnlyDiagonal) 
 TAYdCYsYsCAYd = Matmul2(Transpose(AYd),CYsYsCAYd,OnlyDiagonal) 
 TAYdCYzTYzCAYd = Matmul2(Transpose(AYd),CYzTYzCAYd,OnlyDiagonal) 
 TAYdCAYdTYdCYd = Matmul2(Transpose(AYd),CAYdTYdCYd,OnlyDiagonal) 
 TAYdCAYsYsCYd = Matmul2(Transpose(AYd),CAYsYsCYd,OnlyDiagonal) 
 TAYdCAYzTYzCYd = Matmul2(Transpose(AYd),CAYzTYzCYd,OnlyDiagonal) 
 TAYeCYeTYeCAYe = Matmul2(Transpose(AYe),CYeTYeCAYe,OnlyDiagonal) 
 TAYeCAYeTYeCYe = Matmul2(Transpose(AYe),CAYeTYeCYe,OnlyDiagonal) 
 TAYuCYuTYuCAYu = Matmul2(Transpose(AYu),CYuTYuCAYu,OnlyDiagonal) 
 TAYuCAYuTYuCYu = Matmul2(Transpose(AYu),CAYuTYuCYu,OnlyDiagonal) 
 TAYzCYdTYdCAYz = Matmul2(Transpose(AYz),CYdTYdCAYz,OnlyDiagonal) 
 TAYzCYsYsCAYz = Matmul2(Transpose(AYz),CYsYsCAYz,OnlyDiagonal) 
 TAYzCYzTYzCAYz = Matmul2(Transpose(AYz),CYzTYzCAYz,OnlyDiagonal) 
 TAYzCAYdTYdCYz = Matmul2(Transpose(AYz),CAYdTYdCYz,OnlyDiagonal) 
 TAYzCAYsYsCYz = Matmul2(Transpose(AYz),CAYsYsCYz,OnlyDiagonal) 
 TAYzCAYzTYzCYz = Matmul2(Transpose(AYz),CAYzTYzCYz,OnlyDiagonal) 
 md2YdadjYdYdadjYd = Matmul2(md2,YdadjYdYdadjYd,OnlyDiagonal) 
 md2YdadjYdYsCYs = Matmul2(md2,YdadjYdYsCYs,OnlyDiagonal) 
 md2YdadjYuYuadjYd = Matmul2(md2,YdadjYuYuadjYd,OnlyDiagonal) 
 md2YsCYdTYdCYs = Matmul2(md2,YsCYdTYdCYs,OnlyDiagonal) 
 md2YsCYsYdadjYd = Matmul2(md2,YsCYsYdadjYd,OnlyDiagonal) 
 md2YsCYsYsCYs = Matmul2(md2,YsCYsYsCYs,OnlyDiagonal) 
 md2YsCYsYzadjYz = Matmul2(md2,YsCYsYzadjYz,OnlyDiagonal) 
 md2YsCYzTYzCYs = Matmul2(md2,YsCYzTYzCYs,OnlyDiagonal) 
 md2YzadjYeYeadjYz = Matmul2(md2,YzadjYeYeadjYz,OnlyDiagonal) 
 md2YzadjYzYsCYs = Matmul2(md2,YzadjYzYsCYs,OnlyDiagonal) 
 md2YzadjYzYzadjYz = Matmul2(md2,YzadjYzYzadjYz,OnlyDiagonal) 
 md2YzCYtYtadjYz = Matmul2(md2,YzCYtYtadjYz,OnlyDiagonal) 
 md2CYsYsCYsYs = Matmul2(md2,CYsYsCYsYs,OnlyDiagonal) 
 me2YeadjYeYeadjYe = Matmul2(me2,YeadjYeYeadjYe,OnlyDiagonal) 
 me2YeadjYzYzadjYe = Matmul2(me2,YeadjYzYzadjYe,OnlyDiagonal) 
 me2YeCYtYtadjYe = Matmul2(me2,YeCYtYtadjYe,OnlyDiagonal) 
 ml2YtadjYeYeCYt = Matmul2(ml2,YtadjYeYeCYt,OnlyDiagonal) 
 ml2YtadjYzYzCYt = Matmul2(ml2,YtadjYzYzCYt,OnlyDiagonal) 
 ml2YtCYtYtCYt = Matmul2(ml2,YtCYtYtCYt,OnlyDiagonal) 
 ml2adjYzYsCYsYz = Matmul2(ml2,adjYzYsCYsYz,OnlyDiagonal) 
 ml2CYtYtCYtYt = Matmul2(ml2,CYtYtCYtYt,OnlyDiagonal) 
 ml2TYeCYeTYeCYe = Matmul2(ml2,TYeCYeTYeCYe,OnlyDiagonal) 
 ml2TYzCYdTYdCYz = Matmul2(ml2,TYzCYdTYdCYz,OnlyDiagonal) 
 ml2TYzCYsYsCYz = Matmul2(ml2,TYzCYsYsCYz,OnlyDiagonal) 
 ml2TYzCYzTYzCYz = Matmul2(ml2,TYzCYzTYzCYz,OnlyDiagonal) 
 mq2adjYdYsCYsYd = Matmul2(mq2,adjYdYsCYsYd,OnlyDiagonal) 
 mq2TYdCYdTYdCYd = Matmul2(mq2,TYdCYdTYdCYd,OnlyDiagonal) 
 mq2TYdCYsYsCYd = Matmul2(mq2,TYdCYsYsCYd,OnlyDiagonal) 
 mq2TYdCYzTYzCYd = Matmul2(mq2,TYdCYzTYzCYd,OnlyDiagonal) 
 mq2TYuCYuTYuCYu = Matmul2(mq2,TYuCYuTYuCYu,OnlyDiagonal) 
 mu2YuadjYdYdadjYu = Matmul2(mu2,YuadjYdYdadjYu,OnlyDiagonal) 
 mu2YuadjYuYuadjYu = Matmul2(mu2,YuadjYuYuadjYu,OnlyDiagonal) 
 Ydmq2adjYdYdadjYd = Matmul2(Yd,mq2adjYdYdadjYd,OnlyDiagonal) 
 Ydmq2adjYdYsCYs = Matmul2(Yd,mq2adjYdYsCYs,OnlyDiagonal) 
 Ydmq2adjYdYzadjYz = Matmul2(Yd,mq2adjYdYzadjYz,OnlyDiagonal) 
 Ydmq2adjYuYuadjYd = Matmul2(Yd,mq2adjYuYuadjYd,OnlyDiagonal) 
 YdadjYdmd2YdadjYd = Matmul2(Yd,adjYdmd2YdadjYd,OnlyDiagonal) 
 YdadjYdmd2YsCYs = Matmul2(Yd,adjYdmd2YsCYs,OnlyDiagonal) 
 YdadjYdmd2YzadjYz = Matmul2(Yd,adjYdmd2YzadjYz,OnlyDiagonal) 
 YdadjYdYdmq2adjYd = Matmul2(Yd,adjYdYdmq2adjYd,OnlyDiagonal) 
 YdadjYdYdadjYdmd2 = Matmul2(Yd,adjYdYdadjYdmd2,OnlyDiagonal) 
 YdadjYdYdadjYdYd = Matmul2(Yd,adjYdYdadjYdYd,OnlyDiagonal) 
 YdadjYdYdadjYdYs = Matmul2(Yd,adjYdYdadjYdYs,OnlyDiagonal) 
 YdadjYdYdadjYdYz = Matmul2(Yd,adjYdYdadjYdYz,OnlyDiagonal) 
 YdadjYdYdadjYdAYd = Matmul2(Yd,adjYdYdadjYdAYd,OnlyDiagonal) 
 YdadjYdYdadjYdAYs = Matmul2(Yd,adjYdYdadjYdAYs,OnlyDiagonal) 
 YdadjYdYdadjYdAYz = Matmul2(Yd,adjYdYdadjYdAYz,OnlyDiagonal) 
 YdadjYdYsCYsYd = Matmul2(Yd,adjYdYsCYsYd,OnlyDiagonal) 
 YdadjYdYsCYsAYd = Matmul2(Yd,adjYdYsCYsAYd,OnlyDiagonal) 
 YdadjYdYzml2adjYz = Matmul2(Yd,adjYdYzml2adjYz,OnlyDiagonal) 
 YdadjYdYzadjYzYd = Matmul2(Yd,adjYdYzadjYzYd,OnlyDiagonal) 
 YdadjYdYzadjYzAYd = Matmul2(Yd,adjYdYzadjYzAYd,OnlyDiagonal) 
 YdadjYdAYdadjYdYd = Matmul2(Yd,adjYdAYdadjYdYd,OnlyDiagonal) 
 YdadjYdAYdadjYdYs = Matmul2(Yd,adjYdAYdadjYdYs,OnlyDiagonal) 
 YdadjYdAYdadjYdYz = Matmul2(Yd,adjYdAYdadjYdYz,OnlyDiagonal) 
 YdadjYdAYsCYsYd = Matmul2(Yd,adjYdAYsCYsYd,OnlyDiagonal) 
 YdadjYdAYzadjYzYd = Matmul2(Yd,adjYdAYzadjYzYd,OnlyDiagonal) 
 YdadjYumu2YuadjYd = Matmul2(Yd,adjYumu2YuadjYd,OnlyDiagonal) 
 YdadjYuYumq2adjYd = Matmul2(Yd,adjYuYumq2adjYd,OnlyDiagonal) 
 YdadjYuYuadjYdmd2 = Matmul2(Yd,adjYuYuadjYdmd2,OnlyDiagonal) 
 YdadjYuYuadjYdYd = Matmul2(Yd,adjYuYuadjYdYd,OnlyDiagonal) 
 YdadjYuYuadjYdYs = Matmul2(Yd,adjYuYuadjYdYs,OnlyDiagonal) 
 YdadjYuYuadjYdYz = Matmul2(Yd,adjYuYuadjYdYz,OnlyDiagonal) 
 YdadjYuYuadjYdAYd = Matmul2(Yd,adjYuYuadjYdAYd,OnlyDiagonal) 
 YdadjYuYuadjYdAYs = Matmul2(Yd,adjYuYuadjYdAYs,OnlyDiagonal) 
 YdadjYuYuadjYdAYz = Matmul2(Yd,adjYuYuadjYdAYz,OnlyDiagonal) 
 YdadjYuYuadjYuYu = Matmul2(Yd,adjYuYuadjYuYu,OnlyDiagonal) 
 YdadjYuYuadjYuAYu = Matmul2(Yd,adjYuYuadjYuAYu,OnlyDiagonal) 
 YdadjYuAYuadjYdYd = Matmul2(Yd,adjYuAYuadjYdYd,OnlyDiagonal) 
 YdadjYuAYuadjYdYs = Matmul2(Yd,adjYuAYuadjYdYs,OnlyDiagonal) 
 YdadjYuAYuadjYdYz = Matmul2(Yd,adjYuAYuadjYdYz,OnlyDiagonal) 
 YdadjYuAYuadjYuYu = Matmul2(Yd,adjYuAYuadjYuYu,OnlyDiagonal) 
 Yeml2adjYeYeadjYe = Matmul2(Ye,ml2adjYeYeadjYe,OnlyDiagonal) 
 Yeml2adjYzYzadjYe = Matmul2(Ye,ml2adjYzYzadjYe,OnlyDiagonal) 
 Yeml2CYtYtadjYe = Matmul2(Ye,ml2CYtYtadjYe,OnlyDiagonal) 
 YeadjYeme2YeadjYe = Matmul2(Ye,adjYeme2YeadjYe,OnlyDiagonal) 
 YeadjYeYeml2adjYe = Matmul2(Ye,adjYeYeml2adjYe,OnlyDiagonal) 
 YeadjYeYeadjYeme2 = Matmul2(Ye,adjYeYeadjYeme2,OnlyDiagonal) 
 YeadjYeYeadjYeYe = Matmul2(Ye,adjYeYeadjYeYe,OnlyDiagonal) 
 YeadjYeYeadjYeAYe = Matmul2(Ye,adjYeYeadjYeAYe,OnlyDiagonal) 
 YeadjYeAYeadjYeYe = Matmul2(Ye,adjYeAYeadjYeYe,OnlyDiagonal) 
 YeadjYzmd2YzadjYe = Matmul2(Ye,adjYzmd2YzadjYe,OnlyDiagonal) 
 YeadjYzYdadjYdYz = Matmul2(Ye,adjYzYdadjYdYz,OnlyDiagonal) 
 YeadjYzYdadjYdAYz = Matmul2(Ye,adjYzYdadjYdAYz,OnlyDiagonal) 
 YeadjYzYsCYsYz = Matmul2(Ye,adjYzYsCYsYz,OnlyDiagonal) 
 YeadjYzYsCYsAYz = Matmul2(Ye,adjYzYsCYsAYz,OnlyDiagonal) 
 YeadjYzYzml2adjYe = Matmul2(Ye,adjYzYzml2adjYe,OnlyDiagonal) 
 YeadjYzYzadjYeme2 = Matmul2(Ye,adjYzYzadjYeme2,OnlyDiagonal) 
 YeadjYzYzadjYeYe = Matmul2(Ye,adjYzYzadjYeYe,OnlyDiagonal) 
 YeadjYzYzadjYeAYe = Matmul2(Ye,adjYzYzadjYeAYe,OnlyDiagonal) 
 YeadjYzYzadjYzYz = Matmul2(Ye,adjYzYzadjYzYz,OnlyDiagonal) 
 YeadjYzYzadjYzAYz = Matmul2(Ye,adjYzYzadjYzAYz,OnlyDiagonal) 
 YeadjYzAYdadjYdYz = Matmul2(Ye,adjYzAYdadjYdYz,OnlyDiagonal) 
 YeadjYzAYsCYsYz = Matmul2(Ye,adjYzAYsCYsYz,OnlyDiagonal) 
 YeadjYzAYzadjYeYe = Matmul2(Ye,adjYzAYzadjYeYe,OnlyDiagonal) 
 YeadjYzAYzadjYzYz = Matmul2(Ye,adjYzAYzadjYzYz,OnlyDiagonal) 
 YeCYtml2YtadjYe = Matmul2(Ye,CYtml2YtadjYe,OnlyDiagonal) 
 YeCYtYtml2adjYe = Matmul2(Ye,CYtYtml2adjYe,OnlyDiagonal) 
 YeCYtYtadjYeme2 = Matmul2(Ye,CYtYtadjYeme2,OnlyDiagonal) 
 YeCYtYtadjYeYe = Matmul2(Ye,CYtYtadjYeYe,OnlyDiagonal) 
 YeCYtYtadjYeAYe = Matmul2(Ye,CYtYtadjYeAYe,OnlyDiagonal) 
 YeCYtYtCYtYt = Matmul2(Ye,CYtYtCYtYt,OnlyDiagonal) 
 YeCYtYtCYtAYt = Matmul2(Ye,CYtYtCYtAYt,OnlyDiagonal) 
 YeCYtAYtadjYeYe = Matmul2(Ye,CYtAYtadjYeYe,OnlyDiagonal) 
 YeCYtAYtCYtYt = Matmul2(Ye,CYtAYtCYtYt,OnlyDiagonal) 
 YeCYtTYeCYeYt = Matmul2(Ye,CYtTYeCYeYt,OnlyDiagonal) 
 YeCYtTYeCYeAYt = Matmul2(Ye,CYtTYeCYeAYt,OnlyDiagonal) 
 YeCYtTYzCYzYt = Matmul2(Ye,CYtTYzCYzYt,OnlyDiagonal) 
 YeCYtTYzCYzAYt = Matmul2(Ye,CYtTYzCYzAYt,OnlyDiagonal) 
 YeCYtTAYeCYeYt = Matmul2(Ye,CYtTAYeCYeYt,OnlyDiagonal) 
 YeCYtTAYzCYzYt = Matmul2(Ye,CYtTAYzCYzYt,OnlyDiagonal) 
 Ysmd2CYdTYdCYs = Matmul2(Ys,md2CYdTYdCYs,OnlyDiagonal) 
 Ysmd2CYsYdadjYd = Matmul2(Ys,md2CYsYdadjYd,OnlyDiagonal) 
 Ysmd2CYsYsCYs = Matmul2(Ys,md2CYsYsCYs,OnlyDiagonal) 
 Ysmd2CYsYzadjYz = Matmul2(Ys,md2CYsYzadjYz,OnlyDiagonal) 
 Ysmd2CYzTYzCYs = Matmul2(Ys,md2CYzTYzCYs,OnlyDiagonal) 
 YsCYdmq2TYdCYs = Matmul2(Ys,CYdmq2TYdCYs,OnlyDiagonal) 
 YsCYdTYdmd2CYs = Matmul2(Ys,CYdTYdmd2CYs,OnlyDiagonal) 
 YsCYdTYdCYdTYd = Matmul2(Ys,CYdTYdCYdTYd,OnlyDiagonal) 
 YsCYdTYdCYdTAYd = Matmul2(Ys,CYdTYdCYdTAYd,OnlyDiagonal) 
 YsCYdTYdCYsmd2 = Matmul2(Ys,CYdTYdCYsmd2,OnlyDiagonal) 
 YsCYdTYdCYsYd = Matmul2(Ys,CYdTYdCYsYd,OnlyDiagonal) 
 YsCYdTYdCYsYs = Matmul2(Ys,CYdTYdCYsYs,OnlyDiagonal) 
 YsCYdTYdCYsYz = Matmul2(Ys,CYdTYdCYsYz,OnlyDiagonal) 
 YsCYdTYdCYsAYd = Matmul2(Ys,CYdTYdCYsAYd,OnlyDiagonal) 
 YsCYdTYdCYsAYs = Matmul2(Ys,CYdTYdCYsAYs,OnlyDiagonal) 
 YsCYdTYdCYsAYz = Matmul2(Ys,CYdTYdCYsAYz,OnlyDiagonal) 
 YsCYdTYuCYuTYd = Matmul2(Ys,CYdTYuCYuTYd,OnlyDiagonal) 
 YsCYdTYuCYuTAYd = Matmul2(Ys,CYdTYuCYuTAYd,OnlyDiagonal) 
 YsCYdTAYdCYdTYd = Matmul2(Ys,CYdTAYdCYdTYd,OnlyDiagonal) 
 YsCYdTAYdCYsYd = Matmul2(Ys,CYdTAYdCYsYd,OnlyDiagonal) 
 YsCYdTAYdCYsYs = Matmul2(Ys,CYdTAYdCYsYs,OnlyDiagonal) 
 YsCYdTAYdCYsYz = Matmul2(Ys,CYdTAYdCYsYz,OnlyDiagonal) 
 YsCYdTAYuCYuTYd = Matmul2(Ys,CYdTAYuCYuTYd,OnlyDiagonal) 
 YsCYsmd2YdadjYd = Matmul2(Ys,CYsmd2YdadjYd,OnlyDiagonal) 
 YsCYsmd2YsCYs = Matmul2(Ys,CYsmd2YsCYs,OnlyDiagonal) 
 YsCYsmd2YzadjYz = Matmul2(Ys,CYsmd2YzadjYz,OnlyDiagonal) 
 YsCYsYdmq2adjYd = Matmul2(Ys,CYsYdmq2adjYd,OnlyDiagonal) 
 YsCYsYdadjYdmd2 = Matmul2(Ys,CYsYdadjYdmd2,OnlyDiagonal) 
 YsCYsYdadjYdYs = Matmul2(Ys,CYsYdadjYdYs,OnlyDiagonal) 
 YsCYsYdadjYdAYs = Matmul2(Ys,CYsYdadjYdAYs,OnlyDiagonal) 
 YsCYsYsmd2CYs = Matmul2(Ys,CYsYsmd2CYs,OnlyDiagonal) 
 YsCYsYsCYsmd2 = Matmul2(Ys,CYsYsCYsmd2,OnlyDiagonal) 
 YsCYsYsCYsYd = Matmul2(Ys,CYsYsCYsYd,OnlyDiagonal) 
 YsCYsYsCYsYs = Matmul2(Ys,CYsYsCYsYs,OnlyDiagonal) 
 YsCYsYsCYsYz = Matmul2(Ys,CYsYsCYsYz,OnlyDiagonal) 
 YsCYsYsCYsAYd = Matmul2(Ys,CYsYsCYsAYd,OnlyDiagonal) 
 YsCYsYsCYsAYs = Matmul2(Ys,CYsYsCYsAYs,OnlyDiagonal) 
 YsCYsYsCYsAYz = Matmul2(Ys,CYsYsCYsAYz,OnlyDiagonal) 
 YsCYsYzml2adjYz = Matmul2(Ys,CYsYzml2adjYz,OnlyDiagonal) 
 YsCYsYzadjYzmd2 = Matmul2(Ys,CYsYzadjYzmd2,OnlyDiagonal) 
 YsCYsYzadjYzYs = Matmul2(Ys,CYsYzadjYzYs,OnlyDiagonal) 
 YsCYsYzadjYzAYs = Matmul2(Ys,CYsYzadjYzAYs,OnlyDiagonal) 
 YsCYsAYdadjYdYs = Matmul2(Ys,CYsAYdadjYdYs,OnlyDiagonal) 
 YsCYsAYsCYsYd = Matmul2(Ys,CYsAYsCYsYd,OnlyDiagonal) 
 YsCYsAYsCYsYs = Matmul2(Ys,CYsAYsCYsYs,OnlyDiagonal) 
 YsCYsAYsCYsYz = Matmul2(Ys,CYsAYsCYsYz,OnlyDiagonal) 
 YsCYsAYzadjYzYs = Matmul2(Ys,CYsAYzadjYzYs,OnlyDiagonal) 
 YsCYzml2TYzCYs = Matmul2(Ys,CYzml2TYzCYs,OnlyDiagonal) 
 YsCYzYtCYtTYz = Matmul2(Ys,CYzYtCYtTYz,OnlyDiagonal) 
 YsCYzYtCYtTAYz = Matmul2(Ys,CYzYtCYtTAYz,OnlyDiagonal) 
 YsCYzAYtCYtTYz = Matmul2(Ys,CYzAYtCYtTYz,OnlyDiagonal) 
 YsCYzTYeCYeTYz = Matmul2(Ys,CYzTYeCYeTYz,OnlyDiagonal) 
 YsCYzTYeCYeTAYz = Matmul2(Ys,CYzTYeCYeTAYz,OnlyDiagonal) 
 YsCYzTYzmd2CYs = Matmul2(Ys,CYzTYzmd2CYs,OnlyDiagonal) 
 YsCYzTYzCYsmd2 = Matmul2(Ys,CYzTYzCYsmd2,OnlyDiagonal) 
 YsCYzTYzCYsYd = Matmul2(Ys,CYzTYzCYsYd,OnlyDiagonal) 
 YsCYzTYzCYsYs = Matmul2(Ys,CYzTYzCYsYs,OnlyDiagonal) 
 YsCYzTYzCYsYz = Matmul2(Ys,CYzTYzCYsYz,OnlyDiagonal) 
 YsCYzTYzCYsAYd = Matmul2(Ys,CYzTYzCYsAYd,OnlyDiagonal) 
 YsCYzTYzCYsAYs = Matmul2(Ys,CYzTYzCYsAYs,OnlyDiagonal) 
 YsCYzTYzCYsAYz = Matmul2(Ys,CYzTYzCYsAYz,OnlyDiagonal) 
 YsCYzTYzCYzTYz = Matmul2(Ys,CYzTYzCYzTYz,OnlyDiagonal) 
 YsCYzTYzCYzTAYz = Matmul2(Ys,CYzTYzCYzTAYz,OnlyDiagonal) 
 YsCYzTAYeCYeTYz = Matmul2(Ys,CYzTAYeCYeTYz,OnlyDiagonal) 
 YsCYzTAYzCYsYd = Matmul2(Ys,CYzTAYzCYsYd,OnlyDiagonal) 
 YsCYzTAYzCYsYs = Matmul2(Ys,CYzTAYzCYsYs,OnlyDiagonal) 
 YsCYzTAYzCYsYz = Matmul2(Ys,CYzTAYzCYsYz,OnlyDiagonal) 
 YsCYzTAYzCYzTYz = Matmul2(Ys,CYzTAYzCYzTYz,OnlyDiagonal) 
 Ytml2adjYeYeCYt = Matmul2(Yt,ml2adjYeYeCYt,OnlyDiagonal) 
 Ytml2adjYzYzCYt = Matmul2(Yt,ml2adjYzYzCYt,OnlyDiagonal) 
 Ytml2CYtYtCYt = Matmul2(Yt,ml2CYtYtCYt,OnlyDiagonal) 
 YtadjYeme2YeCYt = Matmul2(Yt,adjYeme2YeCYt,OnlyDiagonal) 
 YtadjYeYeml2CYt = Matmul2(Yt,adjYeYeml2CYt,OnlyDiagonal) 
 YtadjYeYeadjYeYe = Matmul2(Yt,adjYeYeadjYeYe,OnlyDiagonal) 
 YtadjYeYeadjYeAYe = Matmul2(Yt,adjYeYeadjYeAYe,OnlyDiagonal) 
 YtadjYeYeCYtml2 = Matmul2(Yt,adjYeYeCYtml2,OnlyDiagonal) 
 YtadjYeYeCYtYt = Matmul2(Yt,adjYeYeCYtYt,OnlyDiagonal) 
 YtadjYeYeCYtAYt = Matmul2(Yt,adjYeYeCYtAYt,OnlyDiagonal) 
 YtadjYeAYeadjYeYe = Matmul2(Yt,adjYeAYeadjYeYe,OnlyDiagonal) 
 YtadjYeAYeCYtYt = Matmul2(Yt,adjYeAYeCYtYt,OnlyDiagonal) 
 YtadjYzmd2YzCYt = Matmul2(Yt,adjYzmd2YzCYt,OnlyDiagonal) 
 YtadjYzYdadjYdYz = Matmul2(Yt,adjYzYdadjYdYz,OnlyDiagonal) 
 YtadjYzYdadjYdAYz = Matmul2(Yt,adjYzYdadjYdAYz,OnlyDiagonal) 
 YtadjYzYsCYsYz = Matmul2(Yt,adjYzYsCYsYz,OnlyDiagonal) 
 YtadjYzYsCYsAYz = Matmul2(Yt,adjYzYsCYsAYz,OnlyDiagonal) 
 YtadjYzYzml2CYt = Matmul2(Yt,adjYzYzml2CYt,OnlyDiagonal) 
 YtadjYzYzadjYzYz = Matmul2(Yt,adjYzYzadjYzYz,OnlyDiagonal) 
 YtadjYzYzadjYzAYz = Matmul2(Yt,adjYzYzadjYzAYz,OnlyDiagonal) 
 YtadjYzYzCYtml2 = Matmul2(Yt,adjYzYzCYtml2,OnlyDiagonal) 
 YtadjYzYzCYtYt = Matmul2(Yt,adjYzYzCYtYt,OnlyDiagonal) 
 YtadjYzYzCYtAYt = Matmul2(Yt,adjYzYzCYtAYt,OnlyDiagonal) 
 YtadjYzAYdadjYdYz = Matmul2(Yt,adjYzAYdadjYdYz,OnlyDiagonal) 
 YtadjYzAYsCYsYz = Matmul2(Yt,adjYzAYsCYsYz,OnlyDiagonal) 
 YtadjYzAYzadjYzYz = Matmul2(Yt,adjYzAYzadjYzYz,OnlyDiagonal) 
 YtadjYzAYzCYtYt = Matmul2(Yt,adjYzAYzCYtYt,OnlyDiagonal) 
 YtCYtml2YtCYt = Matmul2(Yt,CYtml2YtCYt,OnlyDiagonal) 
 YtCYtYtml2CYt = Matmul2(Yt,CYtYtml2CYt,OnlyDiagonal) 
 YtCYtYtCYtml2 = Matmul2(Yt,CYtYtCYtml2,OnlyDiagonal) 
 YtCYtYtCYtYt = Matmul2(Yt,CYtYtCYtYt,OnlyDiagonal) 
 YtCYtYtCYtAYt = Matmul2(Yt,CYtYtCYtAYt,OnlyDiagonal) 
 YtCYtAYtCYtYt = Matmul2(Yt,CYtAYtCYtYt,OnlyDiagonal) 
 YtCYtTYeCYeYt = Matmul2(Yt,CYtTYeCYeYt,OnlyDiagonal) 
 YtCYtTYeCYeAYt = Matmul2(Yt,CYtTYeCYeAYt,OnlyDiagonal) 
 YtCYtTYzCYzYt = Matmul2(Yt,CYtTYzCYzYt,OnlyDiagonal) 
 YtCYtTYzCYzAYt = Matmul2(Yt,CYtTYzCYzAYt,OnlyDiagonal) 
 YtCYtTAYeCYeYt = Matmul2(Yt,CYtTAYeCYeYt,OnlyDiagonal) 
 YtCYtTAYzCYzYt = Matmul2(Yt,CYtTAYzCYzYt,OnlyDiagonal) 
 Yumq2adjYdYdadjYu = Matmul2(Yu,mq2adjYdYdadjYu,OnlyDiagonal) 
 Yumq2adjYuYuadjYu = Matmul2(Yu,mq2adjYuYuadjYu,OnlyDiagonal) 
 YuadjYdmd2YdadjYu = Matmul2(Yu,adjYdmd2YdadjYu,OnlyDiagonal) 
 YuadjYdYdmq2adjYu = Matmul2(Yu,adjYdYdmq2adjYu,OnlyDiagonal) 
 YuadjYdYdadjYdYd = Matmul2(Yu,adjYdYdadjYdYd,OnlyDiagonal) 
 YuadjYdYdadjYdAYd = Matmul2(Yu,adjYdYdadjYdAYd,OnlyDiagonal) 
 YuadjYdYdadjYumu2 = Matmul2(Yu,adjYdYdadjYumu2,OnlyDiagonal) 
 YuadjYdYdadjYuYu = Matmul2(Yu,adjYdYdadjYuYu,OnlyDiagonal) 
 YuadjYdYdadjYuAYu = Matmul2(Yu,adjYdYdadjYuAYu,OnlyDiagonal) 
 YuadjYdYsCYsYd = Matmul2(Yu,adjYdYsCYsYd,OnlyDiagonal) 
 YuadjYdYsCYsAYd = Matmul2(Yu,adjYdYsCYsAYd,OnlyDiagonal) 
 YuadjYdYzadjYzYd = Matmul2(Yu,adjYdYzadjYzYd,OnlyDiagonal) 
 YuadjYdYzadjYzAYd = Matmul2(Yu,adjYdYzadjYzAYd,OnlyDiagonal) 
 YuadjYdAYdadjYdYd = Matmul2(Yu,adjYdAYdadjYdYd,OnlyDiagonal) 
 YuadjYdAYdadjYuYu = Matmul2(Yu,adjYdAYdadjYuYu,OnlyDiagonal) 
 YuadjYdAYsCYsYd = Matmul2(Yu,adjYdAYsCYsYd,OnlyDiagonal) 
 YuadjYdAYzadjYzYd = Matmul2(Yu,adjYdAYzadjYzYd,OnlyDiagonal) 
 YuadjYumu2YuadjYu = Matmul2(Yu,adjYumu2YuadjYu,OnlyDiagonal) 
 YuadjYuYumq2adjYu = Matmul2(Yu,adjYuYumq2adjYu,OnlyDiagonal) 
 YuadjYuYuadjYumu2 = Matmul2(Yu,adjYuYuadjYumu2,OnlyDiagonal) 
 YuadjYuYuadjYuYu = Matmul2(Yu,adjYuYuadjYuYu,OnlyDiagonal) 
 YuadjYuYuadjYuAYu = Matmul2(Yu,adjYuYuadjYuAYu,OnlyDiagonal) 
 YuadjYuAYuadjYuYu = Matmul2(Yu,adjYuAYuadjYuYu,OnlyDiagonal) 
 Yzml2adjYeYeadjYz = Matmul2(Yz,ml2adjYeYeadjYz,OnlyDiagonal) 
 Yzml2adjYzYdadjYd = Matmul2(Yz,ml2adjYzYdadjYd,OnlyDiagonal) 
 Yzml2adjYzYsCYs = Matmul2(Yz,ml2adjYzYsCYs,OnlyDiagonal) 
 Yzml2adjYzYzadjYz = Matmul2(Yz,ml2adjYzYzadjYz,OnlyDiagonal) 
 Yzml2CYtYtadjYz = Matmul2(Yz,ml2CYtYtadjYz,OnlyDiagonal) 
 YzadjYeme2YeadjYz = Matmul2(Yz,adjYeme2YeadjYz,OnlyDiagonal) 
 YzadjYeYeml2adjYz = Matmul2(Yz,adjYeYeml2adjYz,OnlyDiagonal) 
 YzadjYeYeadjYeYe = Matmul2(Yz,adjYeYeadjYeYe,OnlyDiagonal) 
 YzadjYeYeadjYeAYe = Matmul2(Yz,adjYeYeadjYeAYe,OnlyDiagonal) 
 YzadjYeYeadjYzmd2 = Matmul2(Yz,adjYeYeadjYzmd2,OnlyDiagonal) 
 YzadjYeYeadjYzYd = Matmul2(Yz,adjYeYeadjYzYd,OnlyDiagonal) 
 YzadjYeYeadjYzYs = Matmul2(Yz,adjYeYeadjYzYs,OnlyDiagonal) 
 YzadjYeYeadjYzYz = Matmul2(Yz,adjYeYeadjYzYz,OnlyDiagonal) 
 YzadjYeYeadjYzAYd = Matmul2(Yz,adjYeYeadjYzAYd,OnlyDiagonal) 
 YzadjYeYeadjYzAYs = Matmul2(Yz,adjYeYeadjYzAYs,OnlyDiagonal) 
 YzadjYeYeadjYzAYz = Matmul2(Yz,adjYeYeadjYzAYz,OnlyDiagonal) 
 YzadjYeAYeadjYeYe = Matmul2(Yz,adjYeAYeadjYeYe,OnlyDiagonal) 
 YzadjYeAYeadjYzYd = Matmul2(Yz,adjYeAYeadjYzYd,OnlyDiagonal) 
 YzadjYeAYeadjYzYs = Matmul2(Yz,adjYeAYeadjYzYs,OnlyDiagonal) 
 YzadjYeAYeadjYzYz = Matmul2(Yz,adjYeAYeadjYzYz,OnlyDiagonal) 
 YzadjYzmd2YdadjYd = Matmul2(Yz,adjYzmd2YdadjYd,OnlyDiagonal) 
 YzadjYzmd2YsCYs = Matmul2(Yz,adjYzmd2YsCYs,OnlyDiagonal) 
 YzadjYzmd2YzadjYz = Matmul2(Yz,adjYzmd2YzadjYz,OnlyDiagonal) 
 YzadjYzYdmq2adjYd = Matmul2(Yz,adjYzYdmq2adjYd,OnlyDiagonal) 
 YzadjYzYdadjYdYz = Matmul2(Yz,adjYzYdadjYdYz,OnlyDiagonal) 
 YzadjYzYdadjYdAYz = Matmul2(Yz,adjYzYdadjYdAYz,OnlyDiagonal) 
 YzadjYzYsCYsYz = Matmul2(Yz,adjYzYsCYsYz,OnlyDiagonal) 
 YzadjYzYsCYsAYz = Matmul2(Yz,adjYzYsCYsAYz,OnlyDiagonal) 
 YzadjYzYzml2adjYz = Matmul2(Yz,adjYzYzml2adjYz,OnlyDiagonal) 
 YzadjYzYzadjYzmd2 = Matmul2(Yz,adjYzYzadjYzmd2,OnlyDiagonal) 
 YzadjYzYzadjYzYd = Matmul2(Yz,adjYzYzadjYzYd,OnlyDiagonal) 
 YzadjYzYzadjYzYs = Matmul2(Yz,adjYzYzadjYzYs,OnlyDiagonal) 
 YzadjYzYzadjYzYz = Matmul2(Yz,adjYzYzadjYzYz,OnlyDiagonal) 
 YzadjYzYzadjYzAYd = Matmul2(Yz,adjYzYzadjYzAYd,OnlyDiagonal) 
 YzadjYzYzadjYzAYs = Matmul2(Yz,adjYzYzadjYzAYs,OnlyDiagonal) 
 YzadjYzYzadjYzAYz = Matmul2(Yz,adjYzYzadjYzAYz,OnlyDiagonal) 
 YzadjYzAYdadjYdYz = Matmul2(Yz,adjYzAYdadjYdYz,OnlyDiagonal) 
 YzadjYzAYsCYsYz = Matmul2(Yz,adjYzAYsCYsYz,OnlyDiagonal) 
 YzadjYzAYzadjYzYd = Matmul2(Yz,adjYzAYzadjYzYd,OnlyDiagonal) 
 YzadjYzAYzadjYzYs = Matmul2(Yz,adjYzAYzadjYzYs,OnlyDiagonal) 
 YzadjYzAYzadjYzYz = Matmul2(Yz,adjYzAYzadjYzYz,OnlyDiagonal) 
 YzCYtml2YtadjYz = Matmul2(Yz,CYtml2YtadjYz,OnlyDiagonal) 
 YzCYtYtml2adjYz = Matmul2(Yz,CYtYtml2adjYz,OnlyDiagonal) 
 YzCYtYtadjYzmd2 = Matmul2(Yz,CYtYtadjYzmd2,OnlyDiagonal) 
 YzCYtYtadjYzYd = Matmul2(Yz,CYtYtadjYzYd,OnlyDiagonal) 
 YzCYtYtadjYzYs = Matmul2(Yz,CYtYtadjYzYs,OnlyDiagonal) 
 YzCYtYtadjYzYz = Matmul2(Yz,CYtYtadjYzYz,OnlyDiagonal) 
 YzCYtYtadjYzAYd = Matmul2(Yz,CYtYtadjYzAYd,OnlyDiagonal) 
 YzCYtYtadjYzAYs = Matmul2(Yz,CYtYtadjYzAYs,OnlyDiagonal) 
 YzCYtYtadjYzAYz = Matmul2(Yz,CYtYtadjYzAYz,OnlyDiagonal) 
 YzCYtYtCYtYt = Matmul2(Yz,CYtYtCYtYt,OnlyDiagonal) 
 YzCYtYtCYtAYt = Matmul2(Yz,CYtYtCYtAYt,OnlyDiagonal) 
 YzCYtAYtadjYzYd = Matmul2(Yz,CYtAYtadjYzYd,OnlyDiagonal) 
 YzCYtAYtadjYzYs = Matmul2(Yz,CYtAYtadjYzYs,OnlyDiagonal) 
 YzCYtAYtadjYzYz = Matmul2(Yz,CYtAYtadjYzYz,OnlyDiagonal) 
 YzCYtAYtCYtYt = Matmul2(Yz,CYtAYtCYtYt,OnlyDiagonal) 
 YzCYtTYeCYeYt = Matmul2(Yz,CYtTYeCYeYt,OnlyDiagonal) 
 YzCYtTYeCYeAYt = Matmul2(Yz,CYtTYeCYeAYt,OnlyDiagonal) 
 YzCYtTYzCYzYt = Matmul2(Yz,CYtTYzCYzYt,OnlyDiagonal) 
 YzCYtTYzCYzAYt = Matmul2(Yz,CYtTYzCYzAYt,OnlyDiagonal) 
 YzCYtTAYeCYeYt = Matmul2(Yz,CYtTAYeCYeYt,OnlyDiagonal) 
 YzCYtTAYzCYzYt = Matmul2(Yz,CYtTAYzCYzYt,OnlyDiagonal) 
 AYdadjYdYdadjYdYd = Matmul2(AYd,adjYdYdadjYdYd,OnlyDiagonal) 
 AYdadjYdYdadjYdYs = Matmul2(AYd,adjYdYdadjYdYs,OnlyDiagonal) 
 AYdadjYdYdadjYdYz = Matmul2(AYd,adjYdYdadjYdYz,OnlyDiagonal) 
 AYdadjYdYsCYsYd = Matmul2(AYd,adjYdYsCYsYd,OnlyDiagonal) 
 AYdadjYdYzadjYzYd = Matmul2(AYd,adjYdYzadjYzYd,OnlyDiagonal) 
 AYdadjYuYuadjYdYd = Matmul2(AYd,adjYuYuadjYdYd,OnlyDiagonal) 
 AYdadjYuYuadjYdYs = Matmul2(AYd,adjYuYuadjYdYs,OnlyDiagonal) 
 AYdadjYuYuadjYdYz = Matmul2(AYd,adjYuYuadjYdYz,OnlyDiagonal) 
 AYdadjYuYuadjYuYu = Matmul2(AYd,adjYuYuadjYuYu,OnlyDiagonal) 
 AYeadjYeYeadjYeYe = Matmul2(AYe,adjYeYeadjYeYe,OnlyDiagonal) 
 AYeadjYzYdadjYdYz = Matmul2(AYe,adjYzYdadjYdYz,OnlyDiagonal) 
 AYeadjYzYsCYsYz = Matmul2(AYe,adjYzYsCYsYz,OnlyDiagonal) 
 AYeadjYzYzadjYeYe = Matmul2(AYe,adjYzYzadjYeYe,OnlyDiagonal) 
 AYeadjYzYzadjYzYz = Matmul2(AYe,adjYzYzadjYzYz,OnlyDiagonal) 
 AYeCYtYtadjYeYe = Matmul2(AYe,CYtYtadjYeYe,OnlyDiagonal) 
 AYeCYtYtCYtYt = Matmul2(AYe,CYtYtCYtYt,OnlyDiagonal) 
 AYeCYtTYeCYeYt = Matmul2(AYe,CYtTYeCYeYt,OnlyDiagonal) 
 AYeCYtTYzCYzYt = Matmul2(AYe,CYtTYzCYzYt,OnlyDiagonal) 
 AYsCYdTYdCYdTYd = Matmul2(AYs,CYdTYdCYdTYd,OnlyDiagonal) 
 AYsCYdTYdCYsYd = Matmul2(AYs,CYdTYdCYsYd,OnlyDiagonal) 
 AYsCYdTYdCYsYs = Matmul2(AYs,CYdTYdCYsYs,OnlyDiagonal) 
 AYsCYdTYdCYsYz = Matmul2(AYs,CYdTYdCYsYz,OnlyDiagonal) 
 AYsCYdTYuCYuTYd = Matmul2(AYs,CYdTYuCYuTYd,OnlyDiagonal) 
 AYsCYsYdadjYdYs = Matmul2(AYs,CYsYdadjYdYs,OnlyDiagonal) 
 AYsCYsYsCYsYd = Matmul2(AYs,CYsYsCYsYd,OnlyDiagonal) 
 AYsCYsYsCYsYs = Matmul2(AYs,CYsYsCYsYs,OnlyDiagonal) 
 AYsCYsYsCYsYz = Matmul2(AYs,CYsYsCYsYz,OnlyDiagonal) 
 AYsCYsYzadjYzYs = Matmul2(AYs,CYsYzadjYzYs,OnlyDiagonal) 
 AYsCYzYtCYtTYz = Matmul2(AYs,CYzYtCYtTYz,OnlyDiagonal) 
 AYsCYzTYeCYeTYz = Matmul2(AYs,CYzTYeCYeTYz,OnlyDiagonal) 
 AYsCYzTYzCYsYd = Matmul2(AYs,CYzTYzCYsYd,OnlyDiagonal) 
 AYsCYzTYzCYsYs = Matmul2(AYs,CYzTYzCYsYs,OnlyDiagonal) 
 AYsCYzTYzCYsYz = Matmul2(AYs,CYzTYzCYsYz,OnlyDiagonal) 
 AYsCYzTYzCYzTYz = Matmul2(AYs,CYzTYzCYzTYz,OnlyDiagonal) 
 AYtadjYeYeadjYeYe = Matmul2(AYt,adjYeYeadjYeYe,OnlyDiagonal) 
 AYtadjYeYeCYtYt = Matmul2(AYt,adjYeYeCYtYt,OnlyDiagonal) 
 AYtadjYzYdadjYdYz = Matmul2(AYt,adjYzYdadjYdYz,OnlyDiagonal) 
 AYtadjYzYsCYsYz = Matmul2(AYt,adjYzYsCYsYz,OnlyDiagonal) 
 AYtadjYzYzadjYzYz = Matmul2(AYt,adjYzYzadjYzYz,OnlyDiagonal) 
 AYtadjYzYzCYtYt = Matmul2(AYt,adjYzYzCYtYt,OnlyDiagonal) 
 AYtCYtYtCYtYt = Matmul2(AYt,CYtYtCYtYt,OnlyDiagonal) 
 AYtCYtTYeCYeYt = Matmul2(AYt,CYtTYeCYeYt,OnlyDiagonal) 
 AYtCYtTYzCYzYt = Matmul2(AYt,CYtTYzCYzYt,OnlyDiagonal) 
 AYuadjYdYdadjYdYd = Matmul2(AYu,adjYdYdadjYdYd,OnlyDiagonal) 
 AYuadjYdYdadjYuYu = Matmul2(AYu,adjYdYdadjYuYu,OnlyDiagonal) 
 AYuadjYdYsCYsYd = Matmul2(AYu,adjYdYsCYsYd,OnlyDiagonal) 
 AYuadjYdYzadjYzYd = Matmul2(AYu,adjYdYzadjYzYd,OnlyDiagonal) 
 AYuadjYuYuadjYuYu = Matmul2(AYu,adjYuYuadjYuYu,OnlyDiagonal) 
 AYzadjYeYeadjYeYe = Matmul2(AYz,adjYeYeadjYeYe,OnlyDiagonal) 
 AYzadjYeYeadjYzYd = Matmul2(AYz,adjYeYeadjYzYd,OnlyDiagonal) 
 AYzadjYeYeadjYzYs = Matmul2(AYz,adjYeYeadjYzYs,OnlyDiagonal) 
 AYzadjYeYeadjYzYz = Matmul2(AYz,adjYeYeadjYzYz,OnlyDiagonal) 
 AYzadjYzYdadjYdYz = Matmul2(AYz,adjYzYdadjYdYz,OnlyDiagonal) 
 AYzadjYzYsCYsYz = Matmul2(AYz,adjYzYsCYsYz,OnlyDiagonal) 
 AYzadjYzYzadjYzYd = Matmul2(AYz,adjYzYzadjYzYd,OnlyDiagonal) 
 AYzadjYzYzadjYzYs = Matmul2(AYz,adjYzYzadjYzYs,OnlyDiagonal) 
 AYzadjYzYzadjYzYz = Matmul2(AYz,adjYzYzadjYzYz,OnlyDiagonal) 
 AYzCYtYtadjYzYd = Matmul2(AYz,CYtYtadjYzYd,OnlyDiagonal) 
 AYzCYtYtadjYzYs = Matmul2(AYz,CYtYtadjYzYs,OnlyDiagonal) 
 AYzCYtYtadjYzYz = Matmul2(AYz,CYtYtadjYzYz,OnlyDiagonal) 
 AYzCYtYtCYtYt = Matmul2(AYz,CYtYtCYtYt,OnlyDiagonal) 
 AYzCYtTYeCYeYt = Matmul2(AYz,CYtTYeCYeYt,OnlyDiagonal) 
 AYzCYtTYzCYzYt = Matmul2(AYz,CYtTYzCYzYt,OnlyDiagonal) 
 CYsmd2YsCYsYs = Matmul2(adjYs,md2YsCYsYs,OnlyDiagonal) 
 CYsYdmq2adjYdYs = Matmul2(adjYs,Ydmq2adjYdYs,OnlyDiagonal) 
 CYsYdadjYdmd2Ys = Matmul2(adjYs,YdadjYdmd2Ys,OnlyDiagonal) 
 CYsYdadjYdYsmd2 = Matmul2(adjYs,YdadjYdYsmd2,OnlyDiagonal) 
 CYsYsmd2CYsYs = Matmul2(adjYs,Ysmd2CYsYs,OnlyDiagonal) 
 CYsYsCYsmd2Ys = Matmul2(adjYs,YsCYsmd2Ys,OnlyDiagonal) 
 CYsYsCYsYsmd2 = Matmul2(adjYs,YsCYsYsmd2,OnlyDiagonal) 
 CYsYzml2adjYzYs = Matmul2(adjYs,Yzml2adjYzYs,OnlyDiagonal) 
 CYsYzadjYzmd2Ys = Matmul2(adjYs,YzadjYzmd2Ys,OnlyDiagonal) 
 CYsYzadjYzYsmd2 = Matmul2(adjYs,YzadjYzYsmd2,OnlyDiagonal) 
 CYtml2YtCYtYt = Matmul2(adjYt,ml2YtCYtYt,OnlyDiagonal) 
 CYtYtml2adjYeYe = Matmul2(adjYt,Ytml2adjYeYe,OnlyDiagonal) 
 CYtYtml2adjYzYz = Matmul2(adjYt,Ytml2adjYzYz,OnlyDiagonal) 
 CYtYtml2CYtYt = Matmul2(adjYt,Ytml2CYtYt,OnlyDiagonal) 
 CYtYtadjYeme2Ye = Matmul2(adjYt,YtadjYeme2Ye,OnlyDiagonal) 
 CYtYtadjYeYeml2 = Matmul2(adjYt,YtadjYeYeml2,OnlyDiagonal) 
 CYtYtadjYzmd2Yz = Matmul2(adjYt,YtadjYzmd2Yz,OnlyDiagonal) 
 CYtYtadjYzYzml2 = Matmul2(adjYt,YtadjYzYzml2,OnlyDiagonal) 
 CYtYtCYtml2Yt = Matmul2(adjYt,YtCYtml2Yt,OnlyDiagonal) 
 CYtYtCYtYtml2 = Matmul2(adjYt,YtCYtYtml2,OnlyDiagonal) 
 TYdmd2CYdTYdCYd = Matmul2(Transpose(Yd),md2CYdTYdCYd,OnlyDiagonal) 
 TYdmd2CYsYsCYd = Matmul2(Transpose(Yd),md2CYsYsCYd,OnlyDiagonal) 
 TYdmd2CYzTYzCYd = Matmul2(Transpose(Yd),md2CYzTYzCYd,OnlyDiagonal) 
 TYdCYdmq2TYdCYd = Matmul2(Transpose(Yd),CYdmq2TYdCYd,OnlyDiagonal) 
 TYdCYdTYdmd2CYd = Matmul2(Transpose(Yd),CYdTYdmd2CYd,OnlyDiagonal) 
 TYdCYdTYdCYdmq2 = Matmul2(Transpose(Yd),CYdTYdCYdmq2,OnlyDiagonal) 
 TYdCYsmd2YsCYd = Matmul2(Transpose(Yd),CYsmd2YsCYd,OnlyDiagonal) 
 TYdCYsYsmd2CYd = Matmul2(Transpose(Yd),CYsYsmd2CYd,OnlyDiagonal) 
 TYdCYsYsCYdmq2 = Matmul2(Transpose(Yd),CYsYsCYdmq2,OnlyDiagonal) 
 TYdCYzml2TYzCYd = Matmul2(Transpose(Yd),CYzml2TYzCYd,OnlyDiagonal) 
 TYdCYzTYzmd2CYd = Matmul2(Transpose(Yd),CYzTYzmd2CYd,OnlyDiagonal) 
 TYdCYzTYzCYdmq2 = Matmul2(Transpose(Yd),CYzTYzCYdmq2,OnlyDiagonal) 
 TYeme2CYeTYeCYe = Matmul2(Transpose(Ye),me2CYeTYeCYe,OnlyDiagonal) 
 TYeCYeml2TYeCYe = Matmul2(Transpose(Ye),CYeml2TYeCYe,OnlyDiagonal) 
 TYeCYeTYeme2CYe = Matmul2(Transpose(Ye),CYeTYeme2CYe,OnlyDiagonal) 
 TYeCYeTYeCYeml2 = Matmul2(Transpose(Ye),CYeTYeCYeml2,OnlyDiagonal) 
 TYeCYeTYeCYeYt = Matmul2(Transpose(Ye),CYeTYeCYeYt,OnlyDiagonal) 
 TYeCYeTYeCYeAYt = Matmul2(Transpose(Ye),CYeTYeCYeAYt,OnlyDiagonal) 
 TYeCYeTAYeCYeYt = Matmul2(Transpose(Ye),CYeTAYeCYeYt,OnlyDiagonal) 
 TYumu2CYuTYuCYu = Matmul2(Transpose(Yu),mu2CYuTYuCYu,OnlyDiagonal) 
 TYuCYumq2TYuCYu = Matmul2(Transpose(Yu),CYumq2TYuCYu,OnlyDiagonal) 
 TYuCYuTYumu2CYu = Matmul2(Transpose(Yu),CYuTYumu2CYu,OnlyDiagonal) 
 TYuCYuTYuCYumq2 = Matmul2(Transpose(Yu),CYuTYuCYumq2,OnlyDiagonal) 
 TYzmd2CYdTYdCYz = Matmul2(Transpose(Yz),md2CYdTYdCYz,OnlyDiagonal) 
 TYzmd2CYsYsCYz = Matmul2(Transpose(Yz),md2CYsYsCYz,OnlyDiagonal) 
 TYzmd2CYzTYzCYz = Matmul2(Transpose(Yz),md2CYzTYzCYz,OnlyDiagonal) 
 TYzCYdmq2TYdCYz = Matmul2(Transpose(Yz),CYdmq2TYdCYz,OnlyDiagonal) 
 TYzCYdTYdmd2CYz = Matmul2(Transpose(Yz),CYdTYdmd2CYz,OnlyDiagonal) 
 TYzCYdTYdCYzml2 = Matmul2(Transpose(Yz),CYdTYdCYzml2,OnlyDiagonal) 
 TYzCYdTYdCYzYt = Matmul2(Transpose(Yz),CYdTYdCYzYt,OnlyDiagonal) 
 TYzCYdTYdCYzAYt = Matmul2(Transpose(Yz),CYdTYdCYzAYt,OnlyDiagonal) 
 TYzCYdTAYdCYzYt = Matmul2(Transpose(Yz),CYdTAYdCYzYt,OnlyDiagonal) 
 TYzCYsmd2YsCYz = Matmul2(Transpose(Yz),CYsmd2YsCYz,OnlyDiagonal) 
 TYzCYsYsmd2CYz = Matmul2(Transpose(Yz),CYsYsmd2CYz,OnlyDiagonal) 
 TYzCYsYsCYzml2 = Matmul2(Transpose(Yz),CYsYsCYzml2,OnlyDiagonal) 
 TYzCYsYsCYzYt = Matmul2(Transpose(Yz),CYsYsCYzYt,OnlyDiagonal) 
 TYzCYsYsCYzAYt = Matmul2(Transpose(Yz),CYsYsCYzAYt,OnlyDiagonal) 
 TYzCYsAYsCYzYt = Matmul2(Transpose(Yz),CYsAYsCYzYt,OnlyDiagonal) 
 TYzCYzml2TYzCYz = Matmul2(Transpose(Yz),CYzml2TYzCYz,OnlyDiagonal) 
 TYzCYzTYzmd2CYz = Matmul2(Transpose(Yz),CYzTYzmd2CYz,OnlyDiagonal) 
 TYzCYzTYzCYzml2 = Matmul2(Transpose(Yz),CYzTYzCYzml2,OnlyDiagonal) 
 TYzCYzTYzCYzYt = Matmul2(Transpose(Yz),CYzTYzCYzYt,OnlyDiagonal) 
 TYzCYzTYzCYzAYt = Matmul2(Transpose(Yz),CYzTYzCYzAYt,OnlyDiagonal) 
 TYzCYzTAYzCYzYt = Matmul2(Transpose(Yz),CYzTAYzCYzYt,OnlyDiagonal) 
 TAYeCYeTYeCYeYt = Matmul2(Transpose(AYe),CYeTYeCYeYt,OnlyDiagonal) 
 TAYzCYdTYdCYzYt = Matmul2(Transpose(AYz),CYdTYdCYzYt,OnlyDiagonal) 
 TAYzCYsYsCYzYt = Matmul2(Transpose(AYz),CYsYsCYzYt,OnlyDiagonal) 
 TAYzCYzTYzCYzYt = Matmul2(Transpose(AYz),CYzTYzCYzYt,OnlyDiagonal) 
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
 
betag11 = 68._dp*g13/5._dp

 
If (TwoLoopRGE) Then 
betag12 = 1502._dp*g1**5/75._dp + 174._dp*g13*g22/5._dp + 184._dp*g13*g32/3._dp -& 
&  21._dp*g13*TrCYsYs/5._dp - 9._dp*g13*TrCYtYt/2._dp - 14._dp*g13*TrYdadjYd/5._dp -& 
&  18._dp*g13*TrYeadjYe/5._dp - 3._dp*g13*TrYsCYs/5._dp - 9._dp*g13*TrYtCYt/10._dp -& 
&  26._dp*g13*TrYuadjYu/5._dp - 14._dp*g13*TrYzadjYz/5._dp - 27._dp*g13*L1*Conjg(L1)& 
& /5._dp - 27._dp*g13*L2*Conjg(L2)/5._dp

 
Dg1 = oo16pi2*( betag11 + oo16pi2 * betag12 ) 

 
Else 
Dg1 = oo16pi2* betag11 
End If 
 
 
!-------------------- 
! g2 
!-------------------- 
 
betag21 = 8._dp*g23

 
If (TwoLoopRGE) Then 
betag22 = 58._dp*g12*g23/5._dp + 94._dp*g2**5 + 40._dp*g23*g32 -              & 
&  35._dp*g23*TrCYtYt/6._dp - 6._dp*g23*TrYdadjYd - 2._dp*g23*TrYeadjYe -          & 
&  7._dp*g23*TrYtCYt/6._dp - 6._dp*g23*TrYuadjYu - 6._dp*g23*TrYzadjYz -           & 
&  7._dp*g23*L1*Conjg(L1) - 7._dp*g23*L2*Conjg(L2)

 
Dg2 = oo16pi2*( betag21 + oo16pi2 * betag22 ) 

 
Else 
Dg2 = oo16pi2* betag21 
End If 
 
 
!-------------------- 
! g3 
!-------------------- 
 
betag31 = 4._dp*g33

 
If (TwoLoopRGE) Then 
betag32 = 23._dp*g12*g33/3._dp + 15._dp*g22*g33 + 400._dp*g3**5/3._dp -       & 
&  63._dp*g33*TrCYsYs/8._dp - 4._dp*g33*TrYdadjYd - 9._dp*g33*TrYsCYs/8._dp -      & 
&  4._dp*g33*TrYuadjYu - 4._dp*g33*TrYzadjYz

 
Dg3 = oo16pi2*( betag31 + oo16pi2 * betag32 ) 

 
Else 
Dg3 = oo16pi2* betag31 
End If 
 
 
!-------------------- 
! Yu 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYu1(i1,i2) = -13._dp*g12*Yu(i1,i2)/15._dp - 3._dp*g22*Yu(i1,i2)               & 
&  - 16._dp*g32*Yu(i1,i2)/3._dp + 3._dp*TrYuadjYu*Yu(i1,i2) + 3._dp*L2*Conjg(L2)       & 
& *Yu(i1,i2) + YuadjYdYd(i1,i2) + 3._dp*YuadjYuYu(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYu2(i1,i2) = 5473._dp*g1**4*Yu(i1,i2)/450._dp + g12*g22*Yu(i1,i2)             & 
&  + 57._dp*g2**4*Yu(i1,i2)/2._dp + 136._dp*g12*g32*Yu(i1,i2)/45._dp +               & 
&  8._dp*g22*g32*Yu(i1,i2) + 320._dp*g3**4*Yu(i1,i2)/9._dp - 3._dp*TrYuadjYdYdadjYu*Yu(i1,i2)& 
&  + 4._dp*g12*TrYuadjYu*Yu(i1,i2)/5._dp + 16._dp*g32*TrYuadjYu*Yu(i1,i2)            & 
&  - 9._dp*TrYuadjYuYuadjYu*Yu(i1,i2) + 18._dp*g12*L2*Conjg(L2)*Yu(i1,i2)              & 
& /5._dp + 12._dp*g22*L2*Conjg(L2)*Yu(i1,i2) - 9._dp*L2*TrYuadjYu*Conjg(L2)            & 
& *Yu(i1,i2) - 12._dp*L2**2*Conjg(L2)**2*Yu(i1,i2) + 2._dp*g12*YuadjYdYd(i1,i2)        & 
& /5._dp - 3._dp*TrYdadjYd*YuadjYdYd(i1,i2) - TrYeadjYe*YuadjYdYd(i1,i2) -               & 
&  3._dp*L1*Conjg(L1)*YuadjYdYd(i1,i2) - 2._dp*YuadjYdYdadjYdYd(i1,i2) - 2._dp*YuadjYdYdadjYuYu(i1,i2)& 
&  - 4._dp*YuadjYdYsCYsYd(i1,i2) - 2._dp*YuadjYdYzadjYzYd(i1,i2) + 2._dp*g12*YuadjYuYu(i1,i2)& 
& /5._dp + 6._dp*g22*YuadjYuYu(i1,i2) - 9._dp*TrYuadjYu*YuadjYuYu(i1,i2)               & 
&  - 9._dp*L2*Conjg(L2)*YuadjYuYu(i1,i2) - 4._dp*YuadjYuYuadjYuYu(i1,i2)
 End Do 
End Do 

 
DYu = oo16pi2*( betaYu1 + oo16pi2 * betaYu2 ) 

 
Else 
DYu = oo16pi2* betaYu1 
End If 
 
 
!-------------------- 
! Yd 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYd1(i1,i2) = -7._dp*g12*Yd(i1,i2)/15._dp - 3._dp*g22*Yd(i1,i2) -              & 
&  16._dp*g32*Yd(i1,i2)/3._dp + 3._dp*TrYdadjYd*Yd(i1,i2) + TrYeadjYe*Yd(i1,i2)        & 
&  + 3._dp*L1*Conjg(L1)*Yd(i1,i2) + 3._dp*YdadjYdYd(i1,i2) + YdadjYuYu(i1,i2)            & 
&  + 4._dp*YsCYsYd(i1,i2) + 2._dp*YzadjYzYd(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYd2(i1,i2) = 581._dp*g1**4*Yd(i1,i2)/90._dp + g12*g22*Yd(i1,i2)               & 
&  + 57._dp*g2**4*Yd(i1,i2)/2._dp + 8._dp*g12*g32*Yd(i1,i2)/9._dp + 8._dp*g22*g32*Yd(i1,i2)& 
&  + 320._dp*g3**4*Yd(i1,i2)/9._dp - 2._dp*g12*TrYdadjYd*Yd(i1,i2)/5._dp +             & 
&  16._dp*g32*TrYdadjYd*Yd(i1,i2) - 9._dp*TrYdadjYdYdadjYd*Yd(i1,i2) - 3._dp*TrYdadjYuYuadjYd*Yd(i1,i2)& 
&  + 6._dp*g12*TrYeadjYe*Yd(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*Yd(i1,i2)            & 
&  - 3._dp*TrYeadjYzYzadjYe*Yd(i1,i2) - 3._dp*TrYeCYtYtadjYe*Yd(i1,i2) - 12._dp*TrYsCYsYdadjYd*Yd(i1,i2)& 
&  - 6._dp*TrYzadjYzYdadjYd*Yd(i1,i2) + 18._dp*g12*L1*Conjg(L1)*Yd(i1,i2)              & 
& /5._dp + 12._dp*g22*L1*Conjg(L1)*Yd(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)              & 
& *Yd(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)*Yd(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)& 
& *Yd(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)*Yd(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)       & 
& **2*Yd(i1,i2) + 4._dp*g12*YdadjYdYd(i1,i2)/5._dp + 6._dp*g22*YdadjYdYd(i1,i2)      & 
&  - 9._dp*TrYdadjYd*YdadjYdYd(i1,i2) - 3._dp*TrYeadjYe*YdadjYdYd(i1,i2) -               & 
&  9._dp*L1*Conjg(L1)*YdadjYdYd(i1,i2) - 4._dp*YdadjYdYdadjYdYd(i1,i2) - 4._dp*YdadjYdYsCYsYd(i1,i2)& 
&  - 2._dp*YdadjYdYzadjYzYd(i1,i2) + 4._dp*g12*YdadjYuYu(i1,i2)/5._dp - 3._dp*TrYuadjYu*YdadjYuYu(i1,i2)& 
&  - 3._dp*L2*Conjg(L2)*YdadjYuYu(i1,i2) - 2._dp*YdadjYuYuadjYdYd(i1,i2) -               & 
&  2._dp*YdadjYuYuadjYuYu(i1,i2) - 8._dp*YsCYdTYdCYsYd(i1,i2) + 32._dp*g12*YsCYsYd(i1,i2)& 
& /15._dp + 80._dp*g32*YsCYsYd(i1,i2)/3._dp - 3._dp*TrCYsYs*YsCYsYd(i1,i2)             & 
&  - TrYsCYs*YsCYsYd(i1,i2) - 16._dp*YsCYsYsCYsYd(i1,i2) - 8._dp*YsCYzTYzCYsYd(i1,i2)    & 
&  - 2._dp*YzadjYeYeadjYzYd(i1,i2) + 2._dp*g12*YzadjYzYd(i1,i2)/5._dp + 6._dp*g22*YzadjYzYd(i1,i2)& 
&  - 2._dp*TrYzadjYz*YzadjYzYd(i1,i2) - 6._dp*YzadjYzYzadjYzYd(i1,i2) - 6._dp*YzCYtYtadjYzYd(i1,i2)
 End Do 
End Do 

 
DYd = oo16pi2*( betaYd1 + oo16pi2 * betaYd2 ) 

 
Else 
DYd = oo16pi2* betaYd1 
End If 
 
 
!-------------------- 
! Ye 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYe1(i1,i2) = -9._dp*g12*Ye(i1,i2)/5._dp - 3._dp*g22*Ye(i1,i2) +               & 
&  3._dp*TrYdadjYd*Ye(i1,i2) + TrYeadjYe*Ye(i1,i2) + 3._dp*L1*Conjg(L1)*Ye(i1,i2)        & 
&  + 3._dp*YeadjYeYe(i1,i2) + 3._dp*YeadjYzYz(i1,i2) + 3._dp*YeCYtYt(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYe2(i1,i2) = 261._dp*g1**4*Ye(i1,i2)/10._dp + 9._dp*g12*g22*Ye(i1,i2)         & 
& /5._dp + 57._dp*g2**4*Ye(i1,i2)/2._dp - 2._dp*g12*TrYdadjYd*Ye(i1,i2)/5._dp +        & 
&  16._dp*g32*TrYdadjYd*Ye(i1,i2) - 9._dp*TrYdadjYdYdadjYd*Ye(i1,i2) - 3._dp*TrYdadjYuYuadjYd*Ye(i1,i2)& 
&  + 6._dp*g12*TrYeadjYe*Ye(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*Ye(i1,i2)            & 
&  - 3._dp*TrYeadjYzYzadjYe*Ye(i1,i2) - 3._dp*TrYeCYtYtadjYe*Ye(i1,i2) - 12._dp*TrYsCYsYdadjYd*Ye(i1,i2)& 
&  - 6._dp*TrYzadjYzYdadjYd*Ye(i1,i2) + 18._dp*g12*L1*Conjg(L1)*Ye(i1,i2)              & 
& /5._dp + 12._dp*g22*L1*Conjg(L1)*Ye(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)              & 
& *Ye(i1,i2)/4._dp - 9._dp*L1*TrYdadjYd*Conjg(L1)*Ye(i1,i2) - 3._dp*L1*TrYeadjYe*Conjg(L1)& 
& *Ye(i1,i2) - 3._dp*L1*TrYtCYt*Conjg(L1)*Ye(i1,i2)/4._dp - 12._dp*L1**2*Conjg(L1)       & 
& **2*Ye(i1,i2) + 6._dp*g22*YeadjYeYe(i1,i2) - 9._dp*TrYdadjYd*YeadjYeYe(i1,i2)        & 
&  - 3._dp*TrYeadjYe*YeadjYeYe(i1,i2) - 9._dp*L1*Conjg(L1)*YeadjYeYe(i1,i2)              & 
&  - 4._dp*YeadjYeYeadjYeYe(i1,i2) - 6._dp*YeadjYzYdadjYdYz(i1,i2) - 12._dp*YeadjYzYsCYsYz(i1,i2)& 
&  - 2._dp*g12*YeadjYzYz(i1,i2)/5._dp + 16._dp*g32*YeadjYzYz(i1,i2) - 3._dp*TrYzadjYz*YeadjYzYz(i1,i2)& 
&  - 6._dp*YeadjYzYzadjYeYe(i1,i2) - 6._dp*YeadjYzYzadjYzYz(i1,i2) - 3._dp*YeCYtTYeCYeYt(i1,i2)& 
&  - 9._dp*YeCYtTYzCYzYt(i1,i2) + 18._dp*g12*YeCYtYt(i1,i2)/5._dp + 12._dp*g22*YeCYtYt(i1,i2)& 
&  - 9._dp*TrCYtYt*YeCYtYt(i1,i2)/4._dp - 3._dp*TrYtCYt*YeCYtYt(i1,i2)/4._dp -           & 
&  3._dp*L1*Conjg(L1)*YeCYtYt(i1,i2) - 6._dp*YeCYtYtadjYeYe(i1,i2) - 9._dp*YeCYtYtCYtYt(i1,i2)
 End Do 
End Do 

 
DYe = oo16pi2*( betaYe1 + oo16pi2 * betaYe2 ) 

 
Else 
DYe = oo16pi2* betaYe1 
End If 
 
 
!-------------------- 
! Yt 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYt1(i1,i2) = (sqrt2*TYeCYeYt(i1,i2) + 3._dp*sqrt2*TYzCYzYt(i1,i2)     & 
&  - 9._dp*sqrt2*g12*Yt(i1,i2)/5._dp - 7._dp*sqrt2*g22*Yt(i1,i2)         & 
&  + 3._dp*TrCYtYt*Yt(i1,i2)/(2._dp*sqrt2) + (TrYtCYt*Yt(i1,i2))/(2._dp*sqrt2)& 
&  + sqrt2*L1*Conjg(L1)*Yt(i1,i2) + sqrt2*YtadjYeYe(i1,i2) + 3._dp*sqrt2& 
& *YtadjYzYz(i1,i2) + 6._dp*sqrt2*YtCYtYt(i1,i2))/sqrt2
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYt2(i1,i2) = (-2._dp*sqrt2*TYeCYeTYeCYeYt(i1,i2) + 6._dp*sqrt2        & 
& *g12*TYeCYeYt(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*TYeCYeYt(i1,i2)             & 
&  - sqrt2*TrYeadjYe*TYeCYeYt(i1,i2) - 3._dp*sqrt2*L1*Conjg(L1)              & 
& *TYeCYeYt(i1,i2) - 6._dp*sqrt2*TYzCYdTYdCYzYt(i1,i2) - 12._dp*sqrt2        & 
& *TYzCYsYsCYzYt(i1,i2) - 6._dp*sqrt2*TYzCYzTYzCYzYt(i1,i2) - 2._dp*sqrt2    & 
& *g12*TYzCYzYt(i1,i2)/5._dp + 16._dp*sqrt2*g32*TYzCYzYt(i1,i2) -              & 
&  3._dp*sqrt2*TrYzadjYz*TYzCYzYt(i1,i2) + 417._dp*g1**4*Yt(i1,i2)/(25._dp*sqrt2)& 
&  + 444._dp*sqrt2*g1**4*Yt(i1,i2)/25._dp + 57._dp*sqrt2*g12*g22*Yt(i1,i2)& 
& /5._dp + 57._dp*(g2**4*Yt(i1,i2))/sqrt2 + 48._dp*sqrt2*g2**4*Yt(i1,i2)     & 
&  - 9._dp*g12*TrCYtYt*Yt(i1,i2)/(10._dp*sqrt2) - 3._dp*g22*TrCYtYt*Yt(i1,i2)  & 
& /(2._dp*sqrt2) - sqrt2*TrCYtYtadjYeYe*Yt(i1,i2) - 3._dp*sqrt2        & 
& *TrCYtYtadjYzYz*Yt(i1,i2) - 19._dp*TrCYtYtCYtYt*Yt(i1,i2)/(2._dp*sqrt2)          & 
&  - (TrYeCYtYtadjYe*Yt(i1,i2))/sqrt2 - (TrYtadjYeYeCYt*Yt(i1,i2))/sqrt2     & 
&  - 3._dp*(TrYtadjYzYzCYt*Yt(i1,i2))/sqrt2 - 3._dp*g12*TrYtCYt*Yt(i1,i2)        & 
& /(10._dp*sqrt2) - (g22*TrYtCYt*Yt(i1,i2))/(2._dp*sqrt2) - 5._dp*TrYtCYtYtCYt*Yt(i1,i2)& 
& /(2._dp*sqrt2) - 3._dp*(TrYzCYtYtadjYz*Yt(i1,i2))/sqrt2 - 3._dp*sqrt2& 
& *g12*L1*Conjg(L1)*Yt(i1,i2)/5._dp - sqrt2*g22*L1*Conjg(L1)*Yt(i1,i2)         & 
&  - 6._dp*sqrt2*L1*TrYdadjYd*Conjg(L1)*Yt(i1,i2) - 2._dp*sqrt2              & 
& *L1*TrYeadjYe*Conjg(L1)*Yt(i1,i2) - 6._dp*sqrt2*L1**2*Conjg(L1)**2*Yt(i1,i2)     & 
&  + 6._dp*sqrt2*g12*YtadjYeYe(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*YtadjYeYe(i1,i2)& 
&  - sqrt2*TrYeadjYe*YtadjYeYe(i1,i2) - 3._dp*sqrt2*L1*Conjg(L1)             & 
& *YtadjYeYe(i1,i2) - 2._dp*sqrt2*YtadjYeYeadjYeYe(i1,i2) - 3._dp*sqrt2      & 
& *YtadjYeYeCYtYt(i1,i2) - 6._dp*sqrt2*YtadjYzYdadjYdYz(i1,i2) - 12._dp*sqrt2& 
& *YtadjYzYsCYsYz(i1,i2) - 2._dp*sqrt2*g12*YtadjYzYz(i1,i2)/5._dp +              & 
&  16._dp*sqrt2*g32*YtadjYzYz(i1,i2) - 3._dp*sqrt2*TrYzadjYz*YtadjYzYz(i1,i2)& 
&  - 6._dp*sqrt2*YtadjYzYzadjYzYz(i1,i2) - 9._dp*sqrt2*YtadjYzYzCYtYt(i1,i2) & 
&  - 3._dp*sqrt2*YtCYtTYeCYeYt(i1,i2) - 9._dp*sqrt2*YtCYtTYzCYzYt(i1,i2)     & 
&  + 36._dp*sqrt2*g12*YtCYtYt(i1,i2)/5._dp + 24._dp*sqrt2*g22*YtCYtYt(i1,i2)& 
&  - 9._dp*(TrCYtYt*YtCYtYt(i1,i2))/sqrt2 - 3._dp*(TrYtCYt*YtCYtYt(i1,i2))         & 
& /sqrt2 - 6._dp*sqrt2*L1*Conjg(L1)*YtCYtYt(i1,i2) - 18._dp*sqrt2      & 
& *YtCYtYtCYtYt(i1,i2))/sqrt2
 End Do 
End Do 

 
DYt = oo16pi2*( betaYt1 + oo16pi2 * betaYt2 ) 

 
Else 
DYt = oo16pi2* betaYt1 
End If 
 
 
!-------------------- 
! Ys 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYs1(i1,i2) = (2._dp*sqrt2*YdadjYdYs(i1,i2) - 4._dp*sqrt2              & 
& *g12*Ys(i1,i2)/5._dp - 12._dp*sqrt2*g32*Ys(i1,i2) + 3._dp*TrCYsYs*Ys(i1,i2)  & 
& /(2._dp*sqrt2) + (TrYsCYs*Ys(i1,i2))/(2._dp*sqrt2) + 2._dp*sqrt2     & 
& *YsCYdTYd(i1,i2) + 8._dp*sqrt2*YsCYsYs(i1,i2) + 2._dp*sqrt2*YsCYzTYz(i1,i2)& 
&  + 2._dp*sqrt2*YzadjYzYs(i1,i2))/sqrt2
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYs2(i1,i2) = (-2._dp*sqrt2*YdadjYdYdadjYdYs(i1,i2) + 2._dp*sqrt2      & 
& *g12*YdadjYdYs(i1,i2)/5._dp + 6._dp*sqrt2*g22*YdadjYdYs(i1,i2)               & 
&  - 6._dp*sqrt2*TrYdadjYd*YdadjYdYs(i1,i2) - 2._dp*sqrt2*TrYeadjYe*YdadjYdYs(i1,i2)& 
&  - 6._dp*sqrt2*L1*Conjg(L1)*YdadjYdYs(i1,i2) - 2._dp*sqrt2*YdadjYuYuadjYdYs(i1,i2)& 
&  + 56._dp*sqrt2*g1**4*Ys(i1,i2)/5._dp + 128._dp*sqrt2*g12*g32*Ys(i1,i2)& 
& /15._dp + 320._dp*sqrt2*g3**4*Ys(i1,i2)/3._dp - 2._dp*sqrt2*TrCYsYdadjYdYs*Ys(i1,i2)& 
&  - (sqrt2*g12*TrCYsYs*Ys(i1,i2))/5._dp - sqrt2*g32*TrCYsYs*Ys(i1,i2)   & 
&  - 7._dp*(TrCYsYsCYsYs*Ys(i1,i2))/sqrt2 - 3._dp*sqrt2*TrCYsYsCYsYs*Ys(i1,i2)& 
&  - 2._dp*sqrt2*TrCYsYzadjYzYs*Ys(i1,i2) - sqrt2*TrYdadjYdYsCYs*Ys(i1,i2)   & 
&  - (sqrt2*g12*TrYsCYs*Ys(i1,i2))/15._dp - (sqrt2*g32*TrYsCYs*Ys(i1,i2))& 
& /3._dp - sqrt2*TrYsCYsYdadjYd*Ys(i1,i2) - (TrYsCYsYsCYs*Ys(i1,i2))               & 
& /sqrt2 - sqrt2*TrYsCYsYsCYs*Ys(i1,i2) - sqrt2*TrYsCYsYzadjYz*Ys(i1,i2)& 
&  - sqrt2*TrYzadjYzYsCYs*Ys(i1,i2) + 2._dp*sqrt2*g12*YsCYdTYd(i1,i2)      & 
& /5._dp + 6._dp*sqrt2*g22*YsCYdTYd(i1,i2) - 6._dp*sqrt2*TrYdadjYd*YsCYdTYd(i1,i2)& 
&  - 2._dp*sqrt2*TrYeadjYe*YsCYdTYd(i1,i2) - 6._dp*sqrt2*L1*Conjg(L1)        & 
& *YsCYdTYd(i1,i2) - 2._dp*sqrt2*YsCYdTYdCYdTYd(i1,i2) - 8._dp*sqrt2         & 
& *YsCYdTYdCYsYs(i1,i2) - 2._dp*sqrt2*YsCYdTYuCYuTYd(i1,i2) - 8._dp*sqrt2    & 
& *YsCYsYdadjYdYs(i1,i2) + 64._dp*sqrt2*g12*YsCYsYs(i1,i2)/15._dp +              & 
&  160._dp*sqrt2*g32*YsCYsYs(i1,i2)/3._dp - 6._dp*sqrt2*TrCYsYs*YsCYsYs(i1,i2)& 
&  - 2._dp*sqrt2*TrYsCYs*YsCYsYs(i1,i2) - 32._dp*sqrt2*YsCYsYsCYsYs(i1,i2)   & 
&  - 8._dp*sqrt2*YsCYsYzadjYzYs(i1,i2) - 2._dp*sqrt2*YsCYzTYeCYeTYz(i1,i2)   & 
&  + 2._dp*sqrt2*g12*YsCYzTYz(i1,i2)/5._dp + 6._dp*sqrt2*g22*YsCYzTYz(i1,i2)& 
&  - 2._dp*sqrt2*TrYzadjYz*YsCYzTYz(i1,i2) - 8._dp*sqrt2*YsCYzTYzCYsYs(i1,i2)& 
&  - 6._dp*sqrt2*YsCYzTYzCYzTYz(i1,i2) - 6._dp*sqrt2*YsCYzYtCYtTYz(i1,i2)    & 
&  - 2._dp*sqrt2*YzadjYeYeadjYzYs(i1,i2) + 2._dp*sqrt2*g12*YzadjYzYs(i1,i2)& 
& /5._dp + 6._dp*sqrt2*g22*YzadjYzYs(i1,i2) - 2._dp*sqrt2*TrYzadjYz*YzadjYzYs(i1,i2)& 
&  - 6._dp*sqrt2*YzadjYzYzadjYzYs(i1,i2) - 6._dp*sqrt2*YzCYtYtadjYzYs(i1,i2))& 
& /sqrt2
 End Do 
End Do 

 
DYs = oo16pi2*( betaYs1 + oo16pi2 * betaYs2 ) 

 
Else 
DYs = oo16pi2* betaYs1 
End If 
 
 
!-------------------- 
! Yz 
!-------------------- 
 
Do i1 = 1,3
Do i2 = 1,3
betaYz1(i1,i2) = 2._dp*YdadjYdYz(i1,i2) + 4._dp*YsCYsYz(i1,i2) - 7._dp*g12*Yz(i1,i2)& 
& /15._dp - 3._dp*g22*Yz(i1,i2) - 16._dp*g32*Yz(i1,i2)/3._dp + TrYzadjYz*Yz(i1,i2)   & 
&  + YzadjYeYe(i1,i2) + 5._dp*YzadjYzYz(i1,i2) + 3._dp*YzCYtYt(i1,i2)
 End Do 
End Do 

 
If (TwoLoopRGE) Then 
Do i1 = 1,3
Do i2 = 1,3
betaYz2(i1,i2) = -2._dp*YdadjYdYdadjYdYz(i1,i2) + 2._dp*g12*YdadjYdYz(i1,i2)        & 
& /5._dp + 6._dp*g22*YdadjYdYz(i1,i2) - 6._dp*TrYdadjYd*YdadjYdYz(i1,i2)               & 
&  - 2._dp*TrYeadjYe*YdadjYdYz(i1,i2) - 6._dp*L1*Conjg(L1)*YdadjYdYz(i1,i2)              & 
&  - 2._dp*YdadjYuYuadjYdYz(i1,i2) - 8._dp*YsCYdTYdCYsYz(i1,i2) - 16._dp*YsCYsYsCYsYz(i1,i2)& 
&  + 32._dp*g12*YsCYsYz(i1,i2)/15._dp + 80._dp*g32*YsCYsYz(i1,i2)/3._dp -            & 
&  3._dp*TrCYsYs*YsCYsYz(i1,i2) - TrYsCYs*YsCYsYz(i1,i2) - 8._dp*YsCYzTYzCYsYz(i1,i2)    & 
&  + 581._dp*g1**4*Yz(i1,i2)/90._dp + g12*g22*Yz(i1,i2) + 57._dp*g2**4*Yz(i1,i2)     & 
& /2._dp + 8._dp*g12*g32*Yz(i1,i2)/9._dp + 8._dp*g22*g32*Yz(i1,i2)               & 
&  + 320._dp*g3**4*Yz(i1,i2)/9._dp - 2._dp*TrYdadjYdYzadjYz*Yz(i1,i2) - 4._dp*TrYsCYsYzadjYz*Yz(i1,i2)& 
&  - TrYzadjYeYeadjYz*Yz(i1,i2) + 2._dp*g12*TrYzadjYz*Yz(i1,i2)/5._dp - 5._dp*TrYzadjYzYzadjYz*Yz(i1,i2)& 
&  - 3._dp*TrYzCYtYtadjYz*Yz(i1,i2) + 6._dp*g12*YzadjYeYe(i1,i2)/5._dp -               & 
&  3._dp*TrYdadjYd*YzadjYeYe(i1,i2) - TrYeadjYe*YzadjYeYe(i1,i2) - 3._dp*L1*Conjg(L1)    & 
& *YzadjYeYe(i1,i2) - 2._dp*YzadjYeYeadjYeYe(i1,i2) - 2._dp*YzadjYeYeadjYzYz(i1,i2)      & 
&  - 6._dp*YzadjYzYdadjYdYz(i1,i2) - 12._dp*YzadjYzYsCYsYz(i1,i2) + 6._dp*g22*YzadjYzYz(i1,i2)& 
&  + 16._dp*g32*YzadjYzYz(i1,i2) - 5._dp*TrYzadjYz*YzadjYzYz(i1,i2) - 12._dp*YzadjYzYzadjYzYz(i1,i2)& 
&  - 3._dp*YzCYtTYeCYeYt(i1,i2) - 9._dp*YzCYtTYzCYzYt(i1,i2) + 18._dp*g12*YzCYtYt(i1,i2)& 
& /5._dp + 12._dp*g22*YzCYtYt(i1,i2) - 9._dp*TrCYtYt*YzCYtYt(i1,i2)/4._dp -            & 
&  3._dp*TrYtCYt*YzCYtYt(i1,i2)/4._dp - 3._dp*L1*Conjg(L1)*YzCYtYt(i1,i2) -              & 
&  6._dp*YzCYtYtadjYzYz(i1,i2) - 9._dp*YzCYtYtCYtYt(i1,i2)
 End Do 
End Do 

 
DYz = oo16pi2*( betaYz1 + oo16pi2 * betaYz2 ) 

 
Else 
DYz = oo16pi2* betaYz1 
End If 
 
 
!-------------------- 
! L1 
!-------------------- 
 
betaL11 = (-9._dp*sqrt2*g12*L1/5._dp - 7._dp*sqrt2*g22*L1 + 3._dp*L1*TrCYtYt/(2._dp*sqrt2)& 
&  + 6._dp*sqrt2*L1*TrYdadjYd + 2._dp*sqrt2*L1*TrYeadjYe + (L1*TrYtCYt)      & 
& /(2._dp*sqrt2) + 7._dp*sqrt2*L1**2*Conjg(L1))/sqrt2

 
If (TwoLoopRGE) Then 
betaL12 = (417._dp*g1**4*L1/(25._dp*sqrt2) + 444._dp*sqrt2*g1**4*L1/25._dp + 57._dp*sqrt2& 
& *g12*g22*L1/5._dp + 57._dp*(g2**4*L1)/sqrt2 + 48._dp*sqrt2             & 
& *g2**4*L1 - 9._dp*g12*L1*TrCYtYt/(10._dp*sqrt2) - 3._dp*g22*L1*TrCYtYt/(2._dp*sqrt2)& 
&  - sqrt2*L1*TrCYtYtadjYeYe - 3._dp*sqrt2*L1*TrCYtYtadjYzYz -               & 
&  19._dp*L1*TrCYtYtCYtYt/(2._dp*sqrt2) - 4._dp*sqrt2*g12*L1*TrYdadjYd/5._dp +& 
&  32._dp*sqrt2*g32*L1*TrYdadjYd - 18._dp*sqrt2*L1*TrYdadjYdYdadjYd -      & 
&  6._dp*sqrt2*L1*TrYdadjYuYuadjYd + 12._dp*sqrt2*g12*L1*TrYeadjYe/5._dp - & 
&  6._dp*sqrt2*L1*TrYeadjYeYeadjYe - 6._dp*sqrt2*L1*TrYeadjYzYzadjYe -       & 
&  7._dp*(L1*TrYeCYtYtadjYe)/sqrt2 - 3._dp*sqrt2*L1*TrYeCYtYtadjYe -         & 
&  24._dp*sqrt2*L1*TrYsCYsYdadjYd - (L1*TrYtadjYeYeCYt)/sqrt2 -              & 
&  3._dp*(L1*TrYtadjYzYzCYt)/sqrt2 - 3._dp*g12*L1*TrYtCYt/(10._dp*sqrt2)   & 
&  - (g22*L1*TrYtCYt)/(2._dp*sqrt2) - 5._dp*L1*TrYtCYtYtCYt/(2._dp*sqrt2)  & 
&  - 12._dp*sqrt2*L1*TrYzadjYzYdadjYd - 3._dp*(L1*TrYzCYtYtadjYz)/sqrt2      & 
&  + 33._dp*sqrt2*g12*L1**2*Conjg(L1)/5._dp + 23._dp*sqrt2*g22*L1**2*Conjg(L1)& 
&  - 9._dp*(L1**2*TrCYtYt*Conjg(L1))/sqrt2 - 24._dp*sqrt2*L1**2*TrYdadjYd*Conjg(L1)& 
&  - 8._dp*sqrt2*L1**2*TrYeadjYe*Conjg(L1) - 3._dp*(L1**2*TrYtCYt*Conjg(L1))       & 
& /sqrt2 - 30._dp*sqrt2*L1**3*Conjg(L1)**2)/sqrt2

 
DL1 = oo16pi2*( betaL11 + oo16pi2 * betaL12 ) 

 
Else 
DL1 = oo16pi2* betaL11 
End If 
 
 
!-------------------- 
! L2 
!-------------------- 
 
betaL21 = (-9._dp*sqrt2*g12*L2/5._dp - 7._dp*sqrt2*g22*L2 + 6._dp*sqrt2& 
& *L2*TrYuadjYu + 7._dp*sqrt2*L2**2*Conjg(L2))/sqrt2

 
If (TwoLoopRGE) Then 
betaL22 = (417._dp*g1**4*L2/(25._dp*sqrt2) + 444._dp*sqrt2*g1**4*L2/25._dp + 57._dp*sqrt2& 
& *g12*g22*L2/5._dp + 57._dp*(g2**4*L2)/sqrt2 + 48._dp*sqrt2             & 
& *g2**4*L2 - 6._dp*sqrt2*L2*TrYuadjYdYdadjYu + 8._dp*sqrt2*g12*L2*TrYuadjYu/5._dp +& 
&  32._dp*sqrt2*g32*L2*TrYuadjYu - 18._dp*sqrt2*L2*TrYuadjYuYuadjYu +      & 
&  33._dp*sqrt2*g12*L2**2*Conjg(L2)/5._dp + 23._dp*sqrt2*g22*L2**2*Conjg(L2)& 
&  - 24._dp*sqrt2*L2**2*TrYuadjYu*Conjg(L2) - 30._dp*sqrt2*L2**3*Conjg(L2)   & 
& **2)/sqrt2

 
DL2 = oo16pi2*( betaL21 + oo16pi2 * betaL22 ) 

 
Else 
DL2 = oo16pi2* betaL21 
End If 
 
 
!-------------------- 
! MTM 
!-------------------- 
 
betaMTM1 = -12._dp*g12*MTM/5._dp - 8._dp*g22*MTM + 3._dp*MTM*TrCYtYt/4._dp +      & 
&  (MTM*TrYtCYt)/4._dp + L1*MTM*Conjg(L1) + L2*MTM*Conjg(L2)

 
If (TwoLoopRGE) Then 
betaMTM2 = 888._dp*g1**4*MTM/25._dp + 96._dp*g12*g22*MTM/5._dp + 96._dp*g2**4*MTM -& 
&  9._dp*g12*MTM*TrCYtYt/20._dp - 3._dp*g22*MTM*TrCYtYt/4._dp - MTM*TrCYtYtadjYeYe - & 
&  3._dp*MTM*TrCYtYtadjYzYz - 19._dp*MTM*TrCYtYtCYtYt/4._dp - (MTM*TrYeCYtYtadjYe)       & 
& /2._dp - (MTM*TrYtadjYeYeCYt)/2._dp - 3._dp*MTM*TrYtadjYzYzCYt/2._dp - 3._dp*g12*MTM*TrYtCYt/20._dp -& 
&  (g22*MTM*TrYtCYt)/4._dp - 5._dp*MTM*TrYtCYtYtCYt/4._dp - 3._dp*MTM*TrYzCYtYtadjYz/2._dp -& 
&  3._dp*g12*L1*MTM*Conjg(L1)/5._dp - g22*L1*MTM*Conjg(L1) - 6._dp*L1*MTM*TrYdadjYd*Conjg(L1)& 
&  - 2._dp*L1*MTM*TrYeadjYe*Conjg(L1) - 6._dp*L1**2*MTM*Conjg(L1)**2 - 3._dp*g12*L2*MTM*Conjg(L2)& 
& /5._dp - g22*L2*MTM*Conjg(L2) - 6._dp*L2*MTM*TrYuadjYu*Conjg(L2) - 6._dp*L2**2*MTM*Conjg(L2)**2

 
DMTM = oo16pi2*( betaMTM1 + oo16pi2 * betaMTM2 ) 

 
Else 
DMTM = oo16pi2* betaMTM1 
End If 
 
 
!-------------------- 
! mue 
!-------------------- 
 
betamue1 = -3._dp*g12*mue/5._dp - 3._dp*g22*mue + 3._dp*TrYdadjYd*mue +           & 
&  TrYeadjYe*mue + 3._dp*TrYuadjYu*mue + 3._dp*L1*mue*Conjg(L1) + 3._dp*L2*mue*Conjg(L2)

 
If (TwoLoopRGE) Then 
betamue2 = 417._dp*g1**4*mue/50._dp + 9._dp*g12*g22*mue/5._dp + 57._dp*g2**4*mue/2._dp -& 
&  2._dp*g12*TrYdadjYd*mue/5._dp + 16._dp*g32*TrYdadjYd*mue - 9._dp*TrYdadjYdYdadjYd*mue -& 
&  3._dp*TrYdadjYuYuadjYd*mue + 6._dp*g12*TrYeadjYe*mue/5._dp - 3._dp*TrYeadjYeYeadjYe*mue -& 
&  3._dp*TrYeadjYzYzadjYe*mue - 3._dp*TrYeCYtYtadjYe*mue - 12._dp*TrYsCYsYdadjYd*mue -   & 
&  3._dp*TrYuadjYdYdadjYu*mue + 4._dp*g12*TrYuadjYu*mue/5._dp + 16._dp*g32*TrYuadjYu*mue -& 
&  9._dp*TrYuadjYuYuadjYu*mue - 6._dp*TrYzadjYzYdadjYd*mue + 18._dp*g12*L1*mue*Conjg(L1)& 
& /5._dp + 12._dp*g22*L1*mue*Conjg(L1) - 9._dp*L1*TrCYtYt*mue*Conjg(L1)/4._dp -        & 
&  9._dp*L1*TrYdadjYd*mue*Conjg(L1) - 3._dp*L1*TrYeadjYe*mue*Conjg(L1) - 3._dp*L1*TrYtCYt*mue*Conjg(L1)& 
& /4._dp - 12._dp*L1**2*mue*Conjg(L1)**2 + 18._dp*g12*L2*mue*Conjg(L2)/5._dp +         & 
&  12._dp*g22*L2*mue*Conjg(L2) - 9._dp*L2*TrYuadjYu*mue*Conjg(L2) - 12._dp*L2**2*mue*Conjg(L2)**2

 
Dmue = oo16pi2*( betamue1 + oo16pi2 * betamue2 ) 

 
Else 
Dmue = oo16pi2* betamue1 
End If 
 
 
!-------------------- 
! MZM 
!-------------------- 
 
betaMZM1 = -1._dp*g12*MZM/15._dp - 3._dp*g22*MZM - 16._dp*g32*MZM/3._dp +       & 
&  MZM*TrYzadjYz

 
If (TwoLoopRGE) Then 
betaMZM2 = 409._dp*g1**4*MZM/450._dp + (g12*g22*MZM)/5._dp + 57._dp*g2**4*MZM/2._dp +& 
&  16._dp*g12*g32*MZM/45._dp + 16._dp*g22*g32*MZM + 320._dp*g3**4*MZM/9._dp -    & 
&  2._dp*MZM*TrYdadjYdYzadjYz - 4._dp*MZM*TrYsCYsYzadjYz - MZM*TrYzadjYeYeadjYz +        & 
&  2._dp*g12*MZM*TrYzadjYz/5._dp - 5._dp*MZM*TrYzadjYzYzadjYz - 3._dp*MZM*TrYzCYtYtadjYz

 
DMZM = oo16pi2*( betaMZM1 + oo16pi2 * betaMZM2 ) 

 
Else 
DMZM = oo16pi2* betaMZM1 
End If 
 
 
!-------------------- 
! MSM 
!-------------------- 
 
betaMSM1 = -16._dp*g12*MSM/15._dp - 40._dp*g32*MSM/3._dp + 3._dp*MSM*TrCYsYs/4._dp +& 
&  (MSM*TrYsCYs)/4._dp

 
If (TwoLoopRGE) Then 
betaMSM2 = 3392._dp*g1**4*MSM/225._dp + 128._dp*g12*g32*MSM/9._dp +               & 
&  1280._dp*g3**4*MSM/9._dp - 2._dp*MSM*TrCYsYdadjYdYs - (g12*MSM*TrCYsYs)             & 
& /5._dp - g32*MSM*TrCYsYs - 13._dp*MSM*TrCYsYsCYsYs/2._dp - 2._dp*MSM*TrCYsYzadjYzYs -& 
&  MSM*TrYdadjYdYsCYs - (g12*MSM*TrYsCYs)/15._dp - (g32*MSM*TrYsCYs)/3._dp -         & 
&  MSM*TrYsCYsYdadjYd - 3._dp*MSM*TrYsCYsYsCYs/2._dp - MSM*TrYsCYsYzadjYz -              & 
&  MSM*TrYzadjYzYsCYs

 
DMSM = oo16pi2*( betaMSM1 + oo16pi2 * betaMSM2 ) 

 
Else 
DMSM = oo16pi2* betaMSM1 
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
&  - 16._dp*g32*AYu(i1,i2)/3._dp + 3._dp*TrYuadjYu*AYu(i1,i2) + 3._dp*L2*Conjg(L2)     & 
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
&  - 10946._dp*g1**4*MassB*Yu(i1,i2)/225._dp - 2._dp*g12*g22*MassB*Yu(i1,i2)         & 
&  - 272._dp*g12*g32*MassB*Yu(i1,i2)/45._dp - 272._dp*g12*g32*MassG*Yu(i1,i2)    & 
& /45._dp - 16._dp*g22*g32*MassG*Yu(i1,i2) - 1280._dp*g3**4*MassG*Yu(i1,i2)          & 
& /9._dp - 2._dp*g12*g22*MassWB*Yu(i1,i2) - 114._dp*g2**4*MassWB*Yu(i1,i2)           & 
&  - 16._dp*g22*g32*MassWB*Yu(i1,i2) + 8._dp*g12*TrAYuadjYu*Yu(i1,i2)              & 
& /5._dp + 32._dp*g32*TrAYuadjYu*Yu(i1,i2) - 6._dp*TrYdadjYuAYuadjYd*Yu(i1,i2)         & 
&  - 6._dp*TrYuadjYdAYdadjYu*Yu(i1,i2) - 8._dp*g12*MassB*TrYuadjYu*Yu(i1,i2)           & 
& /5._dp - 32._dp*g32*MassG*TrYuadjYu*Yu(i1,i2) - 36._dp*TrYuadjYuAYuadjYu*Yu(i1,i2)   & 
&  - 36._dp*g12*L2*MassB*Conjg(L2)*Yu(i1,i2)/5._dp - 24._dp*g22*L2*MassWB*Conjg(L2)  & 
& *Yu(i1,i2) - 18._dp*L2*TrAYuadjYu*Conjg(L2)*Yu(i1,i2) + 36._dp*g12*AL2*Conjg(L2)     & 
& *Yu(i1,i2)/5._dp + 24._dp*g22*AL2*Conjg(L2)*Yu(i1,i2) - 18._dp*TrYuadjYu*AL2*Conjg(L2)& 
& *Yu(i1,i2) - 48._dp*L2*AL2*Conjg(L2)**2*Yu(i1,i2) + 4._dp*g12*YuadjYdAYd(i1,i2)      & 
& /5._dp - 6._dp*TrYdadjYd*YuadjYdAYd(i1,i2) - 2._dp*TrYeadjYe*YuadjYdAYd(i1,i2)         & 
&  - 6._dp*L1*Conjg(L1)*YuadjYdAYd(i1,i2) - 4._dp*YuadjYdAYdadjYdYd(i1,i2)               & 
&  - 4._dp*YuadjYdAYdadjYuYu(i1,i2) - 8._dp*YuadjYdAYsCYsYd(i1,i2) - 4._dp*YuadjYdAYzadjYzYd(i1,i2)& 
&  - 4._dp*g12*MassB*YuadjYdYd(i1,i2)/5._dp - 6._dp*TrAYdadjYd*YuadjYdYd(i1,i2)        & 
&  - 2._dp*TrAYeadjYe*YuadjYdYd(i1,i2) - 6._dp*AL1*Conjg(L1)*YuadjYdYd(i1,i2)            & 
&  - 4._dp*YuadjYdYdadjYdAYd(i1,i2) - 2._dp*YuadjYdYdadjYuAYu(i1,i2) - 8._dp*YuadjYdYsCYsAYd(i1,i2)& 
&  - 4._dp*YuadjYdYzadjYzAYd(i1,i2) + 6._dp*g12*YuadjYuAYu(i1,i2)/5._dp +              & 
&  6._dp*g22*YuadjYuAYu(i1,i2) - 12._dp*TrYuadjYu*YuadjYuAYu(i1,i2) - 12._dp*L2*Conjg(L2)& 
& *YuadjYuAYu(i1,i2) - 8._dp*YuadjYuAYuadjYuYu(i1,i2) - 4._dp*g12*MassB*YuadjYuYu(i1,i2)& 
& /5._dp - 12._dp*g22*MassWB*YuadjYuYu(i1,i2) - 18._dp*TrAYuadjYu*YuadjYuYu(i1,i2)     & 
&  - 18._dp*AL2*Conjg(L2)*YuadjYuYu(i1,i2) - 6._dp*YuadjYuYuadjYuAYu(i1,i2)              & 
&  + 5473._dp*g1**4*AYu(i1,i2)/450._dp + g12*g22*AYu(i1,i2) + 57._dp*g2**4*AYu(i1,i2)& 
& /2._dp + 136._dp*g12*g32*AYu(i1,i2)/45._dp + 8._dp*g22*g32*AYu(i1,i2)          & 
&  + 320._dp*g3**4*AYu(i1,i2)/9._dp - 3._dp*TrYdadjYuYuadjYd*AYu(i1,i2) + 4._dp*g12*TrYuadjYu*AYu(i1,i2)& 
& /5._dp + 16._dp*g32*TrYuadjYu*AYu(i1,i2) - 9._dp*TrYuadjYuYuadjYu*AYu(i1,i2)         & 
&  + 18._dp*g12*L2*Conjg(L2)*AYu(i1,i2)/5._dp + 12._dp*g22*L2*Conjg(L2)              & 
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
&  + TrYeadjYe*AYd(i1,i2) + 3._dp*L1*Conjg(L1)*AYd(i1,i2)
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
& *Conjg(L1) - 3._dp*L2*AYdadjYuYu(i1,i2)*Conjg(L2) - 1162._dp*g1**4*MassB*Yd(i1,i2)     & 
& /45._dp - 2._dp*g12*g22*MassB*Yd(i1,i2) - 16._dp*g12*g32*MassB*Yd(i1,i2)       & 
& /9._dp - 16._dp*g12*g32*MassG*Yd(i1,i2)/9._dp - 16._dp*g22*g32*MassG*Yd(i1,i2) & 
&  - 1280._dp*g3**4*MassG*Yd(i1,i2)/9._dp - 2._dp*g12*g22*MassWB*Yd(i1,i2)           & 
&  - 114._dp*g2**4*MassWB*Yd(i1,i2) - 16._dp*g22*g32*MassWB*Yd(i1,i2) -              & 
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
&  - 12._dp*L1*Conjg(L1)*YdadjYdAYd(i1,i2) - 8._dp*YdadjYdAYdadjYdYd(i1,i2)              & 
&  - 8._dp*YdadjYdAYsCYsYd(i1,i2) - 4._dp*YdadjYdAYzadjYzYd(i1,i2) - 8._dp*g12*MassB*YdadjYdYd(i1,i2)& 
& /5._dp - 12._dp*g22*MassWB*YdadjYdYd(i1,i2) - 18._dp*TrAYdadjYd*YdadjYdYd(i1,i2)     & 
&  - 6._dp*TrAYeadjYe*YdadjYdYd(i1,i2) - 18._dp*AL1*Conjg(L1)*YdadjYdYd(i1,i2)           & 
&  - 6._dp*YdadjYdYdadjYdAYd(i1,i2) - 8._dp*YdadjYdYsCYsAYd(i1,i2) - 4._dp*YdadjYdYzadjYzAYd(i1,i2)& 
&  + 8._dp*g12*YdadjYuAYu(i1,i2)/5._dp - 6._dp*TrYuadjYu*YdadjYuAYu(i1,i2)             & 
&  - 6._dp*L2*Conjg(L2)*YdadjYuAYu(i1,i2) - 4._dp*YdadjYuAYuadjYdYd(i1,i2)               & 
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
&  + 581._dp*g1**4*AYd(i1,i2)/90._dp + g12*g22*AYd(i1,i2) + 57._dp*g2**4*AYd(i1,i2)  & 
& /2._dp + 8._dp*g12*g32*AYd(i1,i2)/9._dp + 8._dp*g22*g32*AYd(i1,i2)             & 
&  + 320._dp*g3**4*AYd(i1,i2)/9._dp - 6._dp*TrCYsYdadjYdYs*AYd(i1,i2) - 3._dp*TrCYtYtadjYeYe*AYd(i1,i2)& 
& /2._dp - 2._dp*g12*TrYdadjYd*AYd(i1,i2)/5._dp + 16._dp*g32*TrYdadjYd*AYd(i1,i2)    & 
&  - 9._dp*TrYdadjYdYdadjYd*AYd(i1,i2) - 9._dp*TrYdadjYdYsCYs*AYd(i1,i2)/2._dp -         & 
&  6._dp*TrYdadjYdYzadjYz*AYd(i1,i2) + 6._dp*g12*TrYeadjYe*AYd(i1,i2)/5._dp -          & 
&  3._dp*TrYeadjYeYeadjYe*AYd(i1,i2) - (TrYeCYtYtadjYe*AYd(i1,i2))/2._dp -               & 
&  3._dp*TrYsCYsYdadjYd*AYd(i1,i2)/2._dp - TrYtadjYeYeCYt*AYd(i1,i2) - 3._dp*TrYuadjYdYdadjYu*AYd(i1,i2)& 
&  - 3._dp*TrYzadjYeYeadjYz*AYd(i1,i2) + 18._dp*g12*L1*Conjg(L1)*AYd(i1,i2)            & 
& /5._dp + 12._dp*g22*L1*Conjg(L1)*AYd(i1,i2) - 9._dp*L1*TrCYtYt*Conjg(L1)             & 
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
&  + 3._dp*L1*Conjg(L1)*AYe(i1,i2)
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
&  - 522._dp*g1**4*MassB*Ye(i1,i2)/5._dp - 18._dp*g12*g22*MassB*Ye(i1,i2)            & 
& /5._dp - 18._dp*g12*g22*MassWB*Ye(i1,i2)/5._dp - 114._dp*g2**4*MassWB*Ye(i1,i2)    & 
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
&  - 12._dp*L1*Conjg(L1)*YeadjYeAYe(i1,i2) - 8._dp*YeadjYeAYeadjYeYe(i1,i2)              & 
&  - 12._dp*g22*MassWB*YeadjYeYe(i1,i2) - 18._dp*TrAYdadjYd*YeadjYeYe(i1,i2)           & 
&  - 6._dp*TrAYeadjYe*YeadjYeYe(i1,i2) - 18._dp*AL1*Conjg(L1)*YeadjYeYe(i1,i2)           & 
&  - 6._dp*YeadjYeYeadjYeAYe(i1,i2) - 12._dp*YeadjYzAYdadjYdYz(i1,i2) - 24._dp*YeadjYzAYsCYsYz(i1,i2)& 
&  - 4._dp*g12*YeadjYzAYz(i1,i2)/5._dp + 32._dp*g32*YeadjYzAYz(i1,i2) -              & 
&  6._dp*TrYzadjYz*YeadjYzAYz(i1,i2) - 12._dp*YeadjYzAYzadjYeYe(i1,i2) - 12._dp*YeadjYzAYzadjYzYz(i1,i2)& 
&  - 12._dp*YeadjYzYdadjYdAYz(i1,i2) - 24._dp*YeadjYzYsCYsAYz(i1,i2) + 4._dp*g12*MassB*YeadjYzYz(i1,i2)& 
& /5._dp - 32._dp*g32*MassG*YeadjYzYz(i1,i2) - 6._dp*TrAYzadjYz*YeadjYzYz(i1,i2)       & 
&  - 6._dp*YeadjYzYzadjYeAYe(i1,i2) - 12._dp*YeadjYzYzadjYzAYz(i1,i2) + 36._dp*g12*YeCYtAYt(i1,i2)& 
& /5._dp + 24._dp*g22*YeCYtAYt(i1,i2) - 9._dp*TrCYtYt*YeCYtAYt(i1,i2)/2._dp -          & 
&  3._dp*TrYtCYt*YeCYtAYt(i1,i2)/2._dp - 6._dp*L1*Conjg(L1)*YeCYtAYt(i1,i2)              & 
&  - 12._dp*YeCYtAYtadjYeYe(i1,i2) - 18._dp*YeCYtAYtCYtYt(i1,i2) - 6._dp*YeCYtTAYeCYeYt(i1,i2)& 
&  - 18._dp*YeCYtTAYzCYzYt(i1,i2) - 6._dp*YeCYtTYeCYeAYt(i1,i2) - 18._dp*YeCYtTYzCYzAYt(i1,i2)& 
&  - 36._dp*g12*MassB*YeCYtYt(i1,i2)/5._dp - 24._dp*g22*MassWB*YeCYtYt(i1,i2)        & 
&  - 6._dp*TrCYtAYt*YeCYtYt(i1,i2) - 6._dp*AL1*Conjg(L1)*YeCYtYt(i1,i2) - 6._dp*YeCYtYtadjYeAYe(i1,i2)& 
&  - 18._dp*YeCYtYtCYtAYt(i1,i2) + 261._dp*g1**4*AYe(i1,i2)/10._dp + 9._dp*g12*g22*AYe(i1,i2)& 
& /5._dp + 57._dp*g2**4*AYe(i1,i2)/2._dp - 6._dp*TrCYsYdadjYdYs*AYe(i1,i2)               & 
&  - 3._dp*TrCYtYtadjYeYe*AYe(i1,i2)/2._dp - 2._dp*g12*TrYdadjYd*AYe(i1,i2)            & 
& /5._dp + 16._dp*g32*TrYdadjYd*AYe(i1,i2) - 9._dp*TrYdadjYdYdadjYd*AYe(i1,i2)         & 
&  - 9._dp*TrYdadjYdYsCYs*AYe(i1,i2)/2._dp - 6._dp*TrYdadjYdYzadjYz*AYe(i1,i2)           & 
&  + 6._dp*g12*TrYeadjYe*AYe(i1,i2)/5._dp - 3._dp*TrYeadjYeYeadjYe*AYe(i1,i2)          & 
&  - (TrYeCYtYtadjYe*AYe(i1,i2))/2._dp - 3._dp*TrYsCYsYdadjYd*AYe(i1,i2)/2._dp -         & 
&  TrYtadjYeYeCYt*AYe(i1,i2) - 3._dp*TrYuadjYdYdadjYu*AYe(i1,i2) - 3._dp*TrYzadjYeYeadjYz*AYe(i1,i2)& 
&  + 18._dp*g12*L1*Conjg(L1)*AYe(i1,i2)/5._dp + 12._dp*g22*L1*Conjg(L1)              & 
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
& *L1*Conjg(L1)*AYt(i1,i2))/sqrt2
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
& *TrYeadjYe*TAYeCYeYt(i1,i2) - 6._dp*sqrt2*L1*Conjg(L1)*TAYeCYeYt(i1,i2)          & 
&  - 12._dp*sqrt2*TAYzCYdTYdCYzYt(i1,i2) - 24._dp*sqrt2*TAYzCYsYsCYzYt(i1,i2)& 
&  - 12._dp*sqrt2*TAYzCYzTYzCYzYt(i1,i2) - 4._dp*sqrt2*g12*TAYzCYzYt(i1,i2)& 
& /5._dp + 32._dp*sqrt2*g32*TAYzCYzYt(i1,i2) - 6._dp*sqrt2*TrYzadjYz*TAYzCYzYt(i1,i2)& 
&  + 6._dp*sqrt2*g12*TYeCYeAYt(i1,i2)/5._dp - 3._dp*sqrt2*TrYdadjYd*TYeCYeAYt(i1,i2)& 
&  - sqrt2*TrYeadjYe*TYeCYeAYt(i1,i2) - 3._dp*sqrt2*L1*Conjg(L1)             & 
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
& *TrAYzadjYz*TYzCYzYt(i1,i2) - 522._dp*sqrt2*g1**4*MassB*Yt(i1,i2)/5._dp -        & 
&  114._dp*sqrt2*g12*g22*MassB*Yt(i1,i2)/5._dp - 114._dp*sqrt2           & 
& *g12*g22*MassWB*Yt(i1,i2)/5._dp - 306._dp*sqrt2*g2**4*MassWB*Yt(i1,i2)       & 
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
& *TrYeadjYe*YtadjYeAYe(i1,i2) - 6._dp*sqrt2*L1*Conjg(L1)*YtadjYeAYe(i1,i2)        & 
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
& /(2._dp*sqrt2) - 9._dp*sqrt2*L1*Conjg(L1)*YtCYtAYt(i1,i2) - 36._dp*sqrt2& 
& *YtCYtAYtCYtYt(i1,i2) - 6._dp*sqrt2*YtCYtTAYeCYeYt(i1,i2) - 18._dp*sqrt2   & 
& *YtCYtTAYzCYzYt(i1,i2) - 6._dp*sqrt2*YtCYtTYeCYeAYt(i1,i2) - 18._dp*sqrt2  & 
& *YtCYtTYzCYzAYt(i1,i2) - 72._dp*sqrt2*g12*MassB*YtCYtYt(i1,i2)/5._dp -         & 
&  48._dp*sqrt2*g22*MassWB*YtCYtYt(i1,i2) - 12._dp*sqrt2*TrCYtAYt*YtCYtYt(i1,i2)& 
&  - 12._dp*sqrt2*AL1*Conjg(L1)*YtCYtYt(i1,i2) - 27._dp*sqrt2*YtCYtYtCYtAYt(i1,i2)& 
&  + 261._dp*g1**4*AYt(i1,i2)/(5._dp*sqrt2) + 57._dp*sqrt2*g12*g22*AYt(i1,i2)& 
& /5._dp + 153._dp*(g2**4*AYt(i1,i2))/sqrt2 - 9._dp*g12*TrCYtYt*AYt(i1,i2)       & 
& /(10._dp*sqrt2) - 3._dp*g22*TrCYtYt*AYt(i1,i2)/(2._dp*sqrt2)             & 
&  - 9._dp*(TrCYtYtCYtYt*AYt(i1,i2))/sqrt2 - 2._dp*sqrt2*TrYeCYtYtadjYe*AYt(i1,i2)& 
&  - 3._dp*g12*TrYtCYt*AYt(i1,i2)/(10._dp*sqrt2) - (g22*TrYtCYt*AYt(i1,i2))    & 
& /(2._dp*sqrt2) - 3._dp*(TrYtCYtYtCYt*AYt(i1,i2))/sqrt2 - 6._dp*sqrt2 & 
& *TrYzCYtYtadjYz*AYt(i1,i2) - 3._dp*sqrt2*g12*L1*Conjg(L1)*AYt(i1,i2)           & 
& /5._dp - sqrt2*g22*L1*Conjg(L1)*AYt(i1,i2) - 6._dp*sqrt2*L1*TrYdadjYd*Conjg(L1)& 
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
&  - 6._dp*sqrt2*L1*Conjg(L1)*YdadjYdAYs(i1,i2) - 2._dp*sqrt2*YdadjYdYdadjYdAYs(i1,i2)& 
&  - 4._dp*sqrt2*g12*MassB*YdadjYdYs(i1,i2)/5._dp - 12._dp*sqrt2           & 
& *g22*MassWB*YdadjYdYs(i1,i2) - 12._dp*sqrt2*TrAYdadjYd*YdadjYdYs(i1,i2)        & 
&  - 4._dp*sqrt2*TrAYeadjYe*YdadjYdYs(i1,i2) - 12._dp*sqrt2*AL1*Conjg(L1)    & 
& *YdadjYdYs(i1,i2) - 4._dp*sqrt2*YdadjYuAYuadjYdYs(i1,i2) - 2._dp*sqrt2     & 
& *YdadjYuYuadjYdAYs(i1,i2) - 224._dp*sqrt2*g1**4*MassB*Ys(i1,i2)/5._dp -          & 
&  256._dp*sqrt2*g12*g32*MassB*Ys(i1,i2)/15._dp - 256._dp*sqrt2          & 
& *g12*g32*MassG*Ys(i1,i2)/15._dp - 1280._dp*sqrt2*g3**4*MassG*Ys(i1,i2)       & 
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
&  - 12._dp*sqrt2*L1*Conjg(L1)*YsCYdTAYd(i1,i2) - 4._dp*sqrt2*YsCYdTAYdCYdTYd(i1,i2)& 
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
&  + 56._dp*sqrt2*g1**4*AYs(i1,i2)/5._dp + 128._dp*sqrt2*g12*g32*AYs(i1,i2)& 
& /15._dp + 320._dp*sqrt2*g3**4*AYs(i1,i2)/3._dp - (sqrt2*g12*TrCYsYs*AYs(i1,i2))& 
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
&  - 2._dp*TrYeadjYe*YdadjYdAYz(i1,i2) - 6._dp*L1*Conjg(L1)*YdadjYdAYz(i1,i2)            & 
&  - 2._dp*YdadjYdYdadjYdAYz(i1,i2) - 4._dp*g12*MassB*YdadjYdYz(i1,i2)/5._dp -         & 
&  12._dp*g22*MassWB*YdadjYdYz(i1,i2) - 12._dp*TrAYdadjYd*YdadjYdYz(i1,i2)             & 
&  - 4._dp*TrAYeadjYe*YdadjYdYz(i1,i2) - 12._dp*AL1*Conjg(L1)*YdadjYdYz(i1,i2)           & 
&  - 4._dp*YdadjYuAYuadjYdYz(i1,i2) - 2._dp*YdadjYuYuadjYdAYz(i1,i2) - 16._dp*YsCYdTAYdCYsYz(i1,i2)& 
&  - 8._dp*YsCYdTYdCYsAYz(i1,i2) - 32._dp*YsCYsAYsCYsYz(i1,i2) + 32._dp*g12*YsCYsAYz(i1,i2)& 
& /15._dp + 80._dp*g32*YsCYsAYz(i1,i2)/3._dp - 3._dp*TrCYsYs*YsCYsAYz(i1,i2)           & 
&  - TrYsCYs*YsCYsAYz(i1,i2) - 16._dp*YsCYsYsCYsAYz(i1,i2) - 64._dp*g12*MassB*YsCYsYz(i1,i2)& 
& /15._dp - 160._dp*g32*MassG*YsCYsYz(i1,i2)/3._dp - 8._dp*TrCYsAYs*YsCYsYz(i1,i2)     & 
&  - 16._dp*YsCYzTAYzCYsYz(i1,i2) - 8._dp*YsCYzTYzCYsAYz(i1,i2) - 1162._dp*g1**4*MassB*Yz(i1,i2)& 
& /45._dp - 2._dp*g12*g22*MassB*Yz(i1,i2) - 16._dp*g12*g32*MassB*Yz(i1,i2)       & 
& /9._dp - 16._dp*g12*g32*MassG*Yz(i1,i2)/9._dp - 16._dp*g22*g32*MassG*Yz(i1,i2) & 
&  - 1280._dp*g3**4*MassG*Yz(i1,i2)/9._dp - 2._dp*g12*g22*MassWB*Yz(i1,i2)           & 
&  - 114._dp*g2**4*MassWB*Yz(i1,i2) - 16._dp*g22*g32*MassWB*Yz(i1,i2) +              & 
&  4._dp*g12*TrAYzadjYz*Yz(i1,i2)/5._dp - 4._dp*TrCYsAYzadjYzYs*Yz(i1,i2)              & 
&  - 5._dp*TrCYsYzadjYzAYs*Yz(i1,i2) - 3._dp*TrCYtAYtadjYzYz*Yz(i1,i2) - 4._dp*TrCYtYtadjYzAYz*Yz(i1,i2)& 
&  - 4._dp*TrYdadjYdAYzadjYz*Yz(i1,i2) - 2._dp*TrYeadjYzAYzadjYe*Yz(i1,i2)               & 
&  - 4._dp*TrYsCYsAYzadjYz*Yz(i1,i2) - 2._dp*TrYtadjYzAYzCYt*Yz(i1,i2) - 2._dp*TrYzadjYeAYeadjYz*Yz(i1,i2)& 
&  - 4._dp*g12*MassB*TrYzadjYz*Yz(i1,i2)/5._dp - 4._dp*TrYzadjYzAYdadjYd*Yz(i1,i2)     & 
&  - 3._dp*TrYzadjYzAYsCYs*Yz(i1,i2) - 20._dp*TrYzadjYzAYzadjYz*Yz(i1,i2) -              & 
&  3._dp*TrYzCYtAYtadjYz*Yz(i1,i2) + 12._dp*g12*YzadjYeAYe(i1,i2)/5._dp -              & 
&  6._dp*TrYdadjYd*YzadjYeAYe(i1,i2) - 2._dp*TrYeadjYe*YzadjYeAYe(i1,i2) -               & 
&  6._dp*L1*Conjg(L1)*YzadjYeAYe(i1,i2) - 4._dp*YzadjYeAYeadjYeYe(i1,i2) -               & 
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
& /2._dp - 6._dp*L1*Conjg(L1)*YzCYtAYt(i1,i2) - 12._dp*YzCYtAYtadjYzYz(i1,i2)            & 
&  - 18._dp*YzCYtAYtCYtYt(i1,i2) - 6._dp*YzCYtTAYeCYeYt(i1,i2) - 18._dp*YzCYtTAYzCYzYt(i1,i2)& 
&  - 6._dp*YzCYtTYeCYeAYt(i1,i2) - 18._dp*YzCYtTYzCYzAYt(i1,i2) - 36._dp*g12*MassB*YzCYtYt(i1,i2)& 
& /5._dp - 24._dp*g22*MassWB*YzCYtYt(i1,i2) - 6._dp*TrCYtAYt*YzCYtYt(i1,i2)            & 
&  - 6._dp*AL1*Conjg(L1)*YzCYtYt(i1,i2) - 6._dp*YzCYtYtadjYzAYz(i1,i2) - 18._dp*YzCYtYtCYtAYt(i1,i2)& 
&  + 581._dp*g1**4*AYz(i1,i2)/90._dp + g12*g22*AYz(i1,i2) + 57._dp*g2**4*AYz(i1,i2)  & 
& /2._dp + 8._dp*g12*g32*AYz(i1,i2)/9._dp + 8._dp*g22*g32*AYz(i1,i2)             & 
&  + 320._dp*g3**4*AYz(i1,i2)/9._dp - 2._dp*TrCYsYzadjYzYs*AYz(i1,i2) - 3._dp*TrCYtYtadjYzYz*AYz(i1,i2)& 
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
betaAL12 = (-522._dp*sqrt2*g1**4*L1*MassB/5._dp - 114._dp*sqrt2           & 
& *g12*g22*L1*MassB/5._dp - 114._dp*sqrt2*g12*g22*L1*MassWB/5._dp -        & 
&  306._dp*sqrt2*g2**4*L1*MassWB - 8._dp*sqrt2*g12*L1*TrAYdadjYd/5._dp +   & 
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
&  12._dp*sqrt2*L1*TrYzCYtAYtadjYz + 261._dp*g1**4*AL1/(5._dp*sqrt2)         & 
&  + 57._dp*sqrt2*g12*g22*AL1/5._dp + 153._dp*(g2**4*AL1)/sqrt2          & 
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
betaAL22 = (-522._dp*sqrt2*g1**4*L2*MassB/5._dp - 114._dp*sqrt2           & 
& *g12*g22*L2*MassB/5._dp - 114._dp*sqrt2*g12*g22*L2*MassWB/5._dp -        & 
&  306._dp*sqrt2*g2**4*L2*MassWB + 16._dp*sqrt2*g12*L2*TrAYuadjYu/5._dp +  & 
&  64._dp*sqrt2*g32*L2*TrAYuadjYu - 12._dp*sqrt2*L2*TrYdadjYuAYuadjYd -    & 
&  12._dp*sqrt2*L2*TrYuadjYdAYdadjYu - 16._dp*sqrt2*g12*L2*MassB*TrYuadjYu/5._dp -& 
&  64._dp*sqrt2*g32*L2*MassG*TrYuadjYu - 72._dp*sqrt2*L2*TrYuadjYuAYuadjYu +& 
&  261._dp*g1**4*AL2/(5._dp*sqrt2) + 57._dp*sqrt2*g12*g22*AL2/5._dp +    & 
&  153._dp*(g2**4*AL2)/sqrt2 - 6._dp*sqrt2*TrYdadjYuYuadjYd*AL2 +            & 
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
betaAmue2 = -834._dp*g1**4*MassB*mue/25._dp - 18._dp*g12*g22*MassB*mue/5._dp -    & 
&  18._dp*g12*g22*MassWB*mue/5._dp - 114._dp*g2**4*MassWB*mue - 4._dp*g12*TrAYdadjYd*mue/5._dp +& 
&  32._dp*g32*TrAYdadjYd*mue + 12._dp*g12*TrAYeadjYe*mue/5._dp + 8._dp*g12*TrAYuadjYu*mue/5._dp +& 
&  32._dp*g32*TrAYuadjYu*mue - 12._dp*TrCYsAYdadjYdYs*mue - 15._dp*TrCYsYdadjYdAYs*mue -& 
&  3._dp*TrCYtAYtadjYeYe*mue - 4._dp*TrCYtYtadjYeAYe*mue + 4._dp*g12*MassB*TrYdadjYd*mue/5._dp -& 
&  32._dp*g32*MassG*TrYdadjYd*mue - 36._dp*TrYdadjYdAYdadjYd*mue - 9._dp*TrYdadjYdAYsCYs*mue -& 
&  12._dp*TrYdadjYdAYzadjYz*mue - 12._dp*TrYdadjYuAYuadjYd*mue - 12._dp*g12*MassB*TrYeadjYe*mue/5._dp -& 
&  12._dp*TrYeadjYeAYeadjYe*mue - 6._dp*TrYeadjYzAYzadjYe*mue - 3._dp*TrYeCYtAYtadjYe*mue -& 
&  12._dp*TrYsCYsAYdadjYd*mue - 2._dp*TrYtadjYeAYeCYt*mue - 12._dp*TrYuadjYdAYdadjYu*mue -& 
&  8._dp*g12*MassB*TrYuadjYu*mue/5._dp - 32._dp*g32*MassG*TrYuadjYu*mue -            & 
&  36._dp*TrYuadjYuAYuadjYu*mue - 6._dp*TrYzadjYeAYeadjYz*mue - 12._dp*TrYzadjYzAYdadjYd*mue +& 
&  417._dp*g1**4*Amue/50._dp + 9._dp*g12*g22*Amue/5._dp + 57._dp*g2**4*Amue/2._dp -  & 
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
betaAMTM2 = -3552._dp*g1**4*MassB*MTM/25._dp - 192._dp*g12*g22*MassB*MTM/5._dp -  & 
&  192._dp*g12*g22*MassWB*MTM/5._dp - 384._dp*g2**4*MassWB*MTM - 6._dp*g12*MTM*TrCYtAYt/5._dp -& 
&  2._dp*g22*MTM*TrCYtAYt - 9._dp*MTM*TrCYtAYtCYtYt + 9._dp*g12*MassB*MTM*TrCYtYt/10._dp +& 
&  3._dp*g22*MassWB*MTM*TrCYtYt/2._dp - 4._dp*MTM*TrCYtYtadjYeAYe - 12._dp*MTM*TrCYtYtadjYzAYz -& 
&  12._dp*MTM*TrCYtYtCYtAYt - 4._dp*MTM*TrYeCYtAYtadjYe + 3._dp*g12*MassB*MTM*TrYtCYt/10._dp +& 
&  (g22*MassWB*MTM*TrYtCYt)/2._dp - 3._dp*MTM*TrYtCYtAYtCYt - 12._dp*MTM*TrYzCYtAYtadjYz +& 
&  888._dp*g1**4*AMTM/25._dp + 96._dp*g12*g22*AMTM/5._dp + 96._dp*g2**4*AMTM -       & 
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
betaAMZM2 = -818._dp*g1**4*MassB*MZM/225._dp - 2._dp*g12*g22*MassB*MZM/5._dp -    & 
&  32._dp*g12*g32*MassB*MZM/45._dp - 32._dp*g12*g32*MassG*MZM/45._dp -           & 
&  32._dp*g22*g32*MassG*MZM - 1280._dp*g3**4*MassG*MZM/9._dp - 2._dp*g12*g22*MassWB*MZM/5._dp -& 
&  114._dp*g2**4*MassWB*MZM - 32._dp*g22*g32*MassWB*MZM + 4._dp*g12*MZM*TrAYzadjYz/5._dp -& 
&  4._dp*MZM*TrCYsAYzadjYzYs - 5._dp*MZM*TrCYsYzadjYzAYs - 3._dp*MZM*TrCYtAYtadjYzYz -   & 
&  4._dp*MZM*TrCYtYtadjYzAYz - 4._dp*MZM*TrYdadjYdAYzadjYz - 2._dp*MZM*TrYeadjYzAYzadjYe -& 
&  4._dp*MZM*TrYsCYsAYzadjYz - 2._dp*MZM*TrYtadjYzAYzCYt - 2._dp*MZM*TrYzadjYeAYeadjYz - & 
&  4._dp*g12*MassB*MZM*TrYzadjYz/5._dp - 4._dp*MZM*TrYzadjYzAYdadjYd - 3._dp*MZM*TrYzadjYzAYsCYs -& 
&  20._dp*MZM*TrYzadjYzAYzadjYz - 3._dp*MZM*TrYzCYtAYtadjYz + 409._dp*g1**4*AMZM/450._dp +& 
&  (g12*g22*AMZM)/5._dp + 57._dp*g2**4*AMZM/2._dp + 16._dp*g12*g32*AMZM/45._dp + & 
&  16._dp*g22*g32*AMZM + 320._dp*g3**4*AMZM/9._dp - 2._dp*TrCYsYzadjYzYs*AMZM -      & 
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
betaAMSM2 = -13568._dp*g1**4*MassB*MSM/225._dp - 256._dp*g12*g32*MassB*MSM/9._dp -& 
&  256._dp*g12*g32*MassG*MSM/9._dp - 5120._dp*g3**4*MassG*MSM/9._dp - 8._dp*g12*MSM*TrCYsAYs/15._dp -& 
&  8._dp*g32*MSM*TrCYsAYs/3._dp - 12._dp*MSM*TrCYsAYsCYsYs - 6._dp*MSM*TrCYsYdadjYdAYs +& 
&  2._dp*g12*MassB*MSM*TrCYsYs/5._dp + 2._dp*g32*MassG*MSM*TrCYsYs - 16._dp*MSM*TrCYsYsCYsAYs -& 
&  6._dp*MSM*TrCYsYzadjYzAYs - 2._dp*MSM*TrYdadjYdAYsCYs + 2._dp*g12*MassB*MSM*TrYsCYs/15._dp +& 
&  2._dp*g32*MassG*MSM*TrYsCYs/3._dp - 8._dp*MSM*TrYsCYsAYdadjYd - 4._dp*MSM*TrYsCYsAYsCYs -& 
&  8._dp*MSM*TrYsCYsAYzadjYz - 2._dp*MSM*TrYzadjYzAYsCYs + 3392._dp*g1**4*AMSM/225._dp + & 
&  128._dp*g12*g32*AMSM/9._dp + 1280._dp*g3**4*AMSM/9._dp - (g12*TrCYsYs*AMSM)     & 
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
&  - TrYeadjYe*mq2TYdCYd(i1,i2) - 3._dp*L1*Conjg(L1)*mq2TYdCYd(i1,i2) - 2._dp*mq2TYdCYdTYdCYd(i1,i2)& 
&  - 4._dp*mq2TYdCYsYsCYd(i1,i2) - 2._dp*mq2TYdCYzTYzCYd(i1,i2) + 4._dp*g12*mq2TYuCYu(i1,i2)& 
& /5._dp - 3._dp*TrYuadjYu*mq2TYuCYu(i1,i2) - 3._dp*L2*Conjg(L2)*mq2TYuCYu(i1,i2)        & 
&  - 2._dp*mq2TYuCYuTYuCYu(i1,i2) + 4._dp*g12*TAYdCAYd(i1,i2)/5._dp - 6._dp*TrYdadjYd*TAYdCAYd(i1,i2)& 
&  - 2._dp*TrYeadjYe*TAYdCAYd(i1,i2) - 6._dp*L1*Conjg(L1)*TAYdCAYd(i1,i2) -              & 
&  4._dp*TAYdCAYdTYdCYd(i1,i2) - 8._dp*TAYdCAYsYsCYd(i1,i2) - 4._dp*TAYdCAYzTYzCYd(i1,i2)& 
&  - 6._dp*TrYdadjAYd*TAYdCYd(i1,i2) - 2._dp*TrYeadjAYe*TAYdCYd(i1,i2) - 4._dp*g12*Conjg(MassB)& 
& *TAYdCYd(i1,i2)/5._dp - 6._dp*L1*Conjg(AL1)*TAYdCYd(i1,i2) - 4._dp*TAYdCYdTYdCAYd(i1,i2)& 
&  - 8._dp*TAYdCYsYsCAYd(i1,i2) - 4._dp*TAYdCYzTYzCAYd(i1,i2) + 8._dp*g12*TAYuCAYu(i1,i2)& 
& /5._dp - 6._dp*TrYuadjYu*TAYuCAYu(i1,i2) - 6._dp*L2*Conjg(L2)*TAYuCAYu(i1,i2)          & 
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
&  - TrYeadjYe*TYdCYdmq2(i1,i2) - 3._dp*L1*Conjg(L1)*TYdCYdmq2(i1,i2) - 4._dp*TYdCYdmq2TYdCYd(i1,i2)& 
&  - 4._dp*TYdCYdTAYdCAYd(i1,i2) - 8._dp*mHd2*TYdCYdTYdCYd(i1,i2) - 2._dp*TYdCYdTYdCYdmq2(i1,i2)& 
&  - 4._dp*TYdCYdTYdmd2CYd(i1,i2) - 8._dp*TYdCYsAYsCAYd(i1,i2) - 8._dp*TYdCYsmd2YsCYd(i1,i2)& 
&  - 8._dp*mHd2*TYdCYsYsCYd(i1,i2) - 8._dp*ms2*TYdCYsYsCYd(i1,i2) - 4._dp*TYdCYsYsCYdmq2(i1,i2)& 
&  - 8._dp*TYdCYsYsmd2CYd(i1,i2) - 4._dp*TYdCYzml2TYzCYd(i1,i2) - 4._dp*TYdCYzTAYzCAYd(i1,i2)& 
&  - 4._dp*mHd2*TYdCYzTYzCYd(i1,i2) - 4._dp*mz2*TYdCYzTYzCYd(i1,i2) - 2._dp*TYdCYzTYzCYdmq2(i1,i2)& 
&  - 4._dp*TYdCYzTYzmd2CYd(i1,i2) + 4._dp*g12*TYdmd2CYd(i1,i2)/5._dp - 6._dp*TrYdadjYd*TYdmd2CYd(i1,i2)& 
&  - 2._dp*TrYeadjYe*TYdmd2CYd(i1,i2) - 6._dp*L1*Conjg(L1)*TYdmd2CYd(i1,i2)              & 
&  - 4._dp*TYdmd2CYdTYdCYd(i1,i2) - 8._dp*TYdmd2CYsYsCYd(i1,i2) - 4._dp*TYdmd2CYzTYzCYd(i1,i2)& 
&  - 8._dp*g12*MassB*TYuCAYu(i1,i2)/5._dp - 6._dp*TrAYuadjYu*TYuCAYu(i1,i2)            & 
&  - 6._dp*AL2*Conjg(L2)*TYuCAYu(i1,i2) - 4._dp*TYuCAYuTAYuCYu(i1,i2) + 8._dp*g12*mHu2*TYuCYu(i1,i2)& 
& /5._dp - 6._dp*TrAYuadjAYu*TYuCYu(i1,i2) - 6._dp*Trmu2YuadjYu*TYuCYu(i1,i2)            & 
&  - 12._dp*mHu2*TrYuadjYu*TYuCYu(i1,i2) - 6._dp*TrYumq2adjYu*TYuCYu(i1,i2)              & 
&  - 18._dp*L2*mHu2*Conjg(L2)*TYuCYu(i1,i2) - 6._dp*L2*mtb2*Conjg(L2)*TYuCYu(i1,i2)      & 
&  + 16._dp*g12*MassB*Conjg(MassB)*TYuCYu(i1,i2)/5._dp - 6._dp*AL2*Conjg(AL2)          & 
& *TYuCYu(i1,i2) + 4._dp*g12*TYuCYumq2(i1,i2)/5._dp - 3._dp*TrYuadjYu*TYuCYumq2(i1,i2) & 
&  - 3._dp*L2*Conjg(L2)*TYuCYumq2(i1,i2) - 4._dp*TYuCYumq2TYuCYu(i1,i2) - 4._dp*TYuCYuTAYuCAYu(i1,i2)& 
&  - 8._dp*mHu2*TYuCYuTYuCYu(i1,i2) - 2._dp*TYuCYuTYuCYumq2(i1,i2) - 4._dp*TYuCYuTYumu2CYu(i1,i2)& 
&  + 8._dp*g12*TYumu2CYu(i1,i2)/5._dp - 6._dp*TrYuadjYu*TYumu2CYu(i1,i2)               & 
&  - 6._dp*L2*Conjg(L2)*TYumu2CYu(i1,i2) - 4._dp*TYumu2CYuTYuCYu(i1,i2)

If (i1.eq.i2) Then 
betamq22(i1,i2) = betamq22(i1,i2)+(409._dp*g1**4*MassB*Conjg(MassB)/75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)& 
& /5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)/45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)& 
& /45._dp + (g12*g22*MassWB*Conjg(MassB))/5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)& 
& /45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 32._dp*g22*g32*MassG*Conjg(MassG)& 
&  + 544._dp*g3**4*MassG*Conjg(MassG)/3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG)     & 
&  + (g12*g22*MassB*Conjg(MassWB))/5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB)    & 
&  + 2._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB)   & 
&  + 32._dp*g22*g32*MassWB*Conjg(MassWB) + 2._dp*g1**4*Tr2(1)/15._dp +               & 
&  6._dp*g2**4*Tr2(2) + 32._dp*g3**4*Tr2(3)/3._dp + 4._dp*g12*Tr3(1)/3._dp +           & 
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
&  - TrYeadjYe*ml2TYeCYe(i1,i2) - 3._dp*L1*Conjg(L1)*ml2TYeCYe(i1,i2) - 2._dp*ml2TYeCYeTYeCYe(i1,i2)& 
&  - 6._dp*ml2TYzCYdTYdCYz(i1,i2) - 12._dp*ml2TYzCYsYsCYz(i1,i2) - 2._dp*g12*ml2TYzCYz(i1,i2)& 
& /5._dp + 16._dp*g32*ml2TYzCYz(i1,i2) - 3._dp*TrYzadjYz*ml2TYzCYz(i1,i2)              & 
&  - 6._dp*ml2TYzCYzTYzCYz(i1,i2) - 3._dp*ml2YtadjYeYeCYt(i1,i2) - 9._dp*ml2YtadjYzYzCYt(i1,i2)& 
&  + 18._dp*g12*ml2YtCYt(i1,i2)/5._dp + 12._dp*g22*ml2YtCYt(i1,i2) - 9._dp*TrCYtYt*ml2YtCYt(i1,i2)& 
& /4._dp - 3._dp*TrYtCYt*ml2YtCYt(i1,i2)/4._dp - 3._dp*L1*Conjg(L1)*ml2YtCYt(i1,i2)      & 
&  - 9._dp*ml2YtCYtYtCYt(i1,i2) + 12._dp*g12*TAYeCAYe(i1,i2)/5._dp - 6._dp*TrYdadjYd*TAYeCAYe(i1,i2)& 
&  - 2._dp*TrYeadjYe*TAYeCAYe(i1,i2) - 6._dp*L1*Conjg(L1)*TAYeCAYe(i1,i2) -              & 
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
&  - TrYeadjYe*TYeCYeml2(i1,i2) - 3._dp*L1*Conjg(L1)*TYeCYeml2(i1,i2) - 4._dp*TYeCYeml2TYeCYe(i1,i2)& 
&  - 4._dp*TYeCYeTAYeCAYe(i1,i2) - 8._dp*mHd2*TYeCYeTYeCYe(i1,i2) - 2._dp*TYeCYeTYeCYeml2(i1,i2)& 
&  - 4._dp*TYeCYeTYeme2CYe(i1,i2) + 12._dp*g12*TYeme2CYe(i1,i2)/5._dp - 6._dp*TrYdadjYd*TYeme2CYe(i1,i2)& 
&  - 2._dp*TrYeadjYe*TYeme2CYe(i1,i2) - 6._dp*L1*Conjg(L1)*TYeme2CYe(i1,i2)              & 
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
& /4._dp - 3._dp*TrYtCYt*YtCYtml2(i1,i2)/4._dp - 3._dp*L1*Conjg(L1)*YtCYtml2(i1,i2)      & 
&  - 18._dp*YtCYtml2YtCYt(i1,i2) - 36._dp*mt2*YtCYtYtCYt(i1,i2) - 9._dp*YtCYtYtCYtml2(i1,i2)& 
&  - 18._dp*YtCYtYtml2CYt(i1,i2) - 6._dp*Ytml2adjYeYeCYt(i1,i2) - 18._dp*Ytml2adjYzYzCYt(i1,i2)& 
&  + 36._dp*g12*Ytml2CYt(i1,i2)/5._dp + 24._dp*g22*Ytml2CYt(i1,i2) - 9._dp*TrCYtYt*Ytml2CYt(i1,i2)& 
& /2._dp - 3._dp*TrYtCYt*Ytml2CYt(i1,i2)/2._dp - 6._dp*L1*Conjg(L1)*Ytml2CYt(i1,i2)      & 
&  - 18._dp*Ytml2CYtYtCYt(i1,i2)

If (i1.eq.i2) Then 
betaml22(i1,i2) = betaml22(i1,i2)+(1251._dp*g1**4*MassB*Conjg(MassB)/25._dp + 18._dp*g12*g22*MassB*Conjg(MassB)& 
& /5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)/5._dp + 9._dp*g12*g22*MassB*Conjg(MassWB)& 
& /5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB)& 
&  + 6._dp*g1**4*Tr2(1)/5._dp + 6._dp*g2**4*Tr2(2) - 4._dp*g12*Tr3(1))
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
& **2 - 48._dp*L1**2*mt2*Conjg(L1)**2 + 1251._dp*g1**4*MassB*Conjg(MassB)/25._dp +       & 
&  18._dp*g12*g22*MassB*Conjg(MassB)/5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)   & 
& /5._dp + 4._dp*g12*TrAYdadjYd*Conjg(MassB)/5._dp - 12._dp*g12*TrAYeadjYe*Conjg(MassB)& 
& /5._dp - 8._dp*g12*MassB*TrYdadjYd*Conjg(MassB)/5._dp + 24._dp*g12*MassB*TrYeadjYe*Conjg(MassB)& 
& /5._dp + 72._dp*g12*L1*MassB*Conjg(L1)*Conjg(MassB)/5._dp - 36._dp*g12*AL1*Conjg(L1)& 
& *Conjg(MassB)/5._dp - 32._dp*g32*TrAYdadjYd*Conjg(MassG) + 64._dp*g32*MassG*TrYdadjYd*Conjg(MassG)& 
&  + 9._dp*g12*g22*MassB*Conjg(MassWB)/5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)& 
& /5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB) + 48._dp*g22*L1*MassWB*Conjg(L1)         & 
& *Conjg(MassWB) - 24._dp*g22*AL1*Conjg(L1)*Conjg(MassWB) - 36._dp*g12*L1*MassB*Conjg(AL1)& 
& /5._dp - 24._dp*g22*L1*MassWB*Conjg(AL1) - 18._dp*L1*TrAYdadjYd*Conjg(AL1)           & 
&  - 6._dp*L1*TrAYeadjYe*Conjg(AL1) - 6._dp*L1*TrCYtAYt*Conjg(AL1) + 36._dp*g12*AL1*Conjg(AL1)& 
& /5._dp + 24._dp*g22*AL1*Conjg(AL1) - 9._dp*TrCYtYt*AL1*Conjg(AL1)/2._dp -            & 
&  18._dp*TrYdadjYd*AL1*Conjg(AL1) - 6._dp*TrYeadjYe*AL1*Conjg(AL1) - 3._dp*TrYtCYt*AL1*Conjg(AL1)& 
& /2._dp - 96._dp*L1*AL1*Conjg(L1)*Conjg(AL1) + 6._dp*g1**4*Tr2(1)/5._dp +               & 
&  6._dp*g2**4*Tr2(2) - 4._dp*g12*Tr3(1)

 
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
& **2 + 1251._dp*g1**4*MassB*Conjg(MassB)/25._dp + 18._dp*g12*g22*MassB*Conjg(MassB) & 
& /5._dp + 9._dp*g12*g22*MassWB*Conjg(MassB)/5._dp - 8._dp*g12*TrAYuadjYu*Conjg(MassB)& 
& /5._dp + 16._dp*g12*MassB*TrYuadjYu*Conjg(MassB)/5._dp + 72._dp*g12*L2*MassB*Conjg(L2)& 
& *Conjg(MassB)/5._dp - 36._dp*g12*AL2*Conjg(L2)*Conjg(MassB)/5._dp - 32._dp*g32*TrAYuadjYu*Conjg(MassG)& 
&  + 64._dp*g32*MassG*TrYuadjYu*Conjg(MassG) + 9._dp*g12*g22*MassB*Conjg(MassWB)   & 
& /5._dp + 18._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB)& 
&  + 48._dp*g22*L2*MassWB*Conjg(L2)*Conjg(MassWB) - 24._dp*g22*AL2*Conjg(L2)         & 
& *Conjg(MassWB) - 36._dp*g12*L2*MassB*Conjg(AL2)/5._dp - 24._dp*g22*L2*MassWB*Conjg(AL2)& 
&  - 18._dp*L2*TrAYuadjYu*Conjg(AL2) + 36._dp*g12*AL2*Conjg(AL2)/5._dp +               & 
&  24._dp*g22*AL2*Conjg(AL2) - 18._dp*TrYuadjYu*AL2*Conjg(AL2) - 96._dp*L2*AL2*Conjg(L2)& 
& *Conjg(AL2) + 6._dp*g1**4*Tr2(1)/5._dp + 6._dp*g2**4*Tr2(2) + 4._dp*g12*Tr3(1)

 
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
&  - 6._dp*L1*Conjg(L1)*md2YdadjYd(i1,i2) - 2._dp*md2YdadjYdYdadjYd(i1,i2)               & 
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
&  - 2._dp*TrYeadjYe*YdadjYdmd2(i1,i2) - 6._dp*L1*Conjg(L1)*YdadjYdmd2(i1,i2)            & 
&  - 4._dp*YdadjYdmd2YdadjYd(i1,i2) - 8._dp*mHd2*YdadjYdYdadjYd(i1,i2) - 2._dp*YdadjYdYdadjYdmd2(i1,i2)& 
&  - 4._dp*YdadjYdYdmq2adjYd(i1,i2) - 4._dp*YdadjYuAYuadjAYd(i1,i2) - 4._dp*YdadjYumu2YuadjYd(i1,i2)& 
&  - 4._dp*mHd2*YdadjYuYuadjYd(i1,i2) - 4._dp*mHu2*YdadjYuYuadjYd(i1,i2) -               & 
&  2._dp*YdadjYuYuadjYdmd2(i1,i2) - 4._dp*YdadjYuYumq2adjYd(i1,i2) + 4._dp*g12*Ydmq2adjYd(i1,i2)& 
& /5._dp + 12._dp*g22*Ydmq2adjYd(i1,i2) - 12._dp*TrYdadjYd*Ydmq2adjYd(i1,i2)           & 
&  - 4._dp*TrYeadjYe*Ydmq2adjYd(i1,i2) - 12._dp*L1*Conjg(L1)*Ydmq2adjYd(i1,i2)           & 
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
betamd22(i1,i2) = betamd22(i1,i2)+(1648._dp*g1**4*MassB*Conjg(MassB)/75._dp + 128._dp*g12*g32*MassB*Conjg(MassB)& 
& /45._dp + 64._dp*g12*g32*MassG*Conjg(MassB)/45._dp + 64._dp*g12*g32*MassB*Conjg(MassG)& 
& /45._dp + 128._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 544._dp*g3**4*MassG*Conjg(MassG)& 
& /3._dp + 8._dp*g1**4*Tr2(1)/15._dp + 32._dp*g3**4*Tr2(3)/3._dp + 8._dp*g12*Tr3(1)    & 
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
&  6._dp*TrYuadjYu*mu2YuadjYu(i1,i2) - 6._dp*L2*Conjg(L2)*mu2YuadjYu(i1,i2)              & 
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
&  6._dp*TrYuadjYu*YuadjYumu2(i1,i2) - 6._dp*L2*Conjg(L2)*YuadjYumu2(i1,i2)              & 
&  - 4._dp*YuadjYumu2YuadjYu(i1,i2) - 8._dp*mHu2*YuadjYuYuadjYu(i1,i2) - 2._dp*YuadjYuYuadjYumu2(i1,i2)& 
&  - 4._dp*YuadjYuYumq2adjYu(i1,i2) - 4._dp*Yumq2adjYdYdadjYu(i1,i2) - 4._dp*g12*Yumq2adjYu(i1,i2)& 
& /5._dp + 12._dp*g22*Yumq2adjYu(i1,i2) - 12._dp*TrYuadjYu*Yumq2adjYu(i1,i2)           & 
&  - 12._dp*L2*Conjg(L2)*Yumq2adjYu(i1,i2) - 4._dp*Yumq2adjYuYuadjYu(i1,i2)

If (i1.eq.i2) Then 
betamu22(i1,i2) = betamu22(i1,i2)+(6784._dp*g1**4*MassB*Conjg(MassB)/75._dp + 512._dp*g12*g32*MassB*Conjg(MassB)& 
& /45._dp + 256._dp*g12*g32*MassG*Conjg(MassB)/45._dp + 256._dp*g12*g32*MassB*Conjg(MassG)& 
& /45._dp + 512._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 544._dp*g3**4*MassG*Conjg(MassG)& 
& /3._dp + 32._dp*g1**4*Tr2(1)/15._dp + 32._dp*g3**4*Tr2(3)/3._dp - 16._dp*g12*Tr3(1)  & 
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
&  - 2._dp*TrYeadjYe*me2YeadjYe(i1,i2) - 6._dp*L1*Conjg(L1)*me2YeadjYe(i1,i2)            & 
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
&  6._dp*L1*Conjg(L1)*YeadjYeme2(i1,i2) - 4._dp*YeadjYeme2YeadjYe(i1,i2) -               & 
&  8._dp*mHd2*YeadjYeYeadjYe(i1,i2) - 2._dp*YeadjYeYeadjYeme2(i1,i2) - 4._dp*YeadjYeYeml2adjYe(i1,i2)& 
&  - 12._dp*YeadjYzAYzadjAYe(i1,i2) - 12._dp*YeadjYzmd2YzadjYe(i1,i2) - 12._dp*mHd2*YeadjYzYzadjYe(i1,i2)& 
&  - 12._dp*mz2*YeadjYzYzadjYe(i1,i2) - 6._dp*YeadjYzYzadjYeme2(i1,i2) - 12._dp*YeadjYzYzml2adjYe(i1,i2)& 
&  - 12._dp*YeCAYtAYtadjYe(i1,i2) - 12._dp*YeCYtAYtadjAYe(i1,i2) - 12._dp*YeCYtml2YtadjYe(i1,i2)& 
&  - 12._dp*mHd2*YeCYtYtadjYe(i1,i2) - 12._dp*mt2*YeCYtYtadjYe(i1,i2) - 6._dp*YeCYtYtadjYeme2(i1,i2)& 
&  - 12._dp*YeCYtYtml2adjYe(i1,i2) - 12._dp*g12*Yeml2adjYe(i1,i2)/5._dp +              & 
&  12._dp*g22*Yeml2adjYe(i1,i2) - 12._dp*TrYdadjYd*Yeml2adjYe(i1,i2) - 4._dp*TrYeadjYe*Yeml2adjYe(i1,i2)& 
&  - 12._dp*L1*Conjg(L1)*Yeml2adjYe(i1,i2) - 4._dp*Yeml2adjYeYeadjYe(i1,i2)              & 
&  - 12._dp*Yeml2adjYzYzadjYe(i1,i2) - 12._dp*Yeml2CYtYtadjYe(i1,i2)

If (i1.eq.i2) Then 
betame22(i1,i2) = betame22(i1,i2)+(5328._dp*g1**4*MassB*Conjg(MassB)/25._dp + 24._dp*g1**4*Tr2(1)& 
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
& **2 + 5328._dp*g1**4*MassB*Conjg(MassB)/25._dp + 192._dp*g12*g22*MassB*Conjg(MassB)& 
& /5._dp + 96._dp*g12*g22*MassWB*Conjg(MassB)/5._dp + 6._dp*g12*TrCYtAYt*Conjg(MassB)& 
& /5._dp - 9._dp*g12*MassB*TrCYtYt*Conjg(MassB)/5._dp - 3._dp*g12*MassB*TrYtCYt*Conjg(MassB)& 
& /5._dp - 12._dp*g12*L1*MassB*Conjg(L1)*Conjg(MassB)/5._dp + 6._dp*g12*AL1*Conjg(L1)& 
& *Conjg(MassB)/5._dp + 96._dp*g12*g22*MassB*Conjg(MassWB)/5._dp + 192._dp*g12*g22*MassWB*Conjg(MassWB)& 
& /5._dp + 544._dp*g2**4*MassWB*Conjg(MassWB) + 2._dp*g22*TrCYtAYt*Conjg(MassWB)       & 
&  - 3._dp*g22*MassWB*TrCYtYt*Conjg(MassWB) - g22*MassWB*TrYtCYt*Conjg(MassWB)       & 
&  - 4._dp*g22*L1*MassWB*Conjg(L1)*Conjg(MassWB) + 2._dp*g22*AL1*Conjg(L1)           & 
& *Conjg(MassWB) + 6._dp*g12*L1*MassB*Conjg(AL1)/5._dp + 2._dp*g22*L1*MassWB*Conjg(AL1)& 
&  - 12._dp*L1*TrAYdadjYd*Conjg(AL1) - 4._dp*L1*TrAYeadjYe*Conjg(AL1) - 6._dp*g12*AL1*Conjg(AL1)& 
& /5._dp - 2._dp*g22*AL1*Conjg(AL1) - 12._dp*TrYdadjYd*AL1*Conjg(AL1) - 4._dp*TrYeadjYe*AL1*Conjg(AL1)& 
&  - 48._dp*L1*AL1*Conjg(L1)*Conjg(AL1) + 24._dp*g1**4*Tr2(1)/5._dp + 16._dp*g2**4*Tr2(2)& 
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
& **2 + 5328._dp*g1**4*MassB*Conjg(MassB)/25._dp + 192._dp*g12*g22*MassB*Conjg(MassB)& 
& /5._dp + 96._dp*g12*g22*MassWB*Conjg(MassB)/5._dp - 12._dp*g12*L2*MassB*Conjg(L2)& 
& *Conjg(MassB)/5._dp + 6._dp*g12*AL2*Conjg(L2)*Conjg(MassB)/5._dp + 96._dp*g12*g22*MassB*Conjg(MassWB)& 
& /5._dp + 192._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 544._dp*g2**4*MassWB*Conjg(MassWB)& 
&  - 4._dp*g22*L2*MassWB*Conjg(L2)*Conjg(MassWB) + 2._dp*g22*AL2*Conjg(L2)           & 
& *Conjg(MassWB) + 6._dp*g12*L2*MassB*Conjg(AL2)/5._dp + 2._dp*g22*L2*MassWB*Conjg(AL2)& 
&  - 12._dp*L2*TrAYuadjYu*Conjg(AL2) - 6._dp*g12*AL2*Conjg(AL2)/5._dp - 2._dp*g22*AL2*Conjg(AL2)& 
&  - 12._dp*TrYuadjYu*AL2*Conjg(AL2) - 48._dp*L2*AL2*Conjg(L2)*Conjg(AL2) +              & 
&  24._dp*g1**4*Tr2(1)/5._dp + 16._dp*g2**4*Tr2(2) - 8._dp*g12*Tr3(1)

 
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
&  2._dp*TrYzadjAYzAYsCYs - 8._dp*TrYzadjYzAYsCAYs - ms2*TrYzadjYzYsCYs + 6784._dp*g1**4*MassB*Conjg(MassB)& 
& /75._dp + 256._dp*g12*g32*MassB*Conjg(MassB)/9._dp + 128._dp*g12*g32*MassG*Conjg(MassB)& 
& /9._dp + 8._dp*g12*TrCYsAYs*Conjg(MassB)/15._dp - 4._dp*g12*MassB*TrCYsYs*Conjg(MassB)& 
& /5._dp - 4._dp*g12*MassB*TrYsCYs*Conjg(MassB)/15._dp + 128._dp*g12*g32*MassB*Conjg(MassG)& 
& /9._dp + 256._dp*g12*g32*MassG*Conjg(MassG)/9._dp + 2320._dp*g3**4*MassG*Conjg(MassG)& 
& /3._dp + 8._dp*g32*TrCYsAYs*Conjg(MassG)/3._dp - 4._dp*g32*MassG*TrCYsYs*Conjg(MassG)& 
&  - 4._dp*g32*MassG*TrYsCYs*Conjg(MassG)/3._dp + 32._dp*g1**4*Tr2(1)/15._dp +         & 
&  80._dp*g3**4*Tr2(3)/3._dp - 16._dp*g12*Tr3(1)/3._dp + 4._dp*g32*Tr3(3)

 
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
betamsb22 = 6784._dp*g1**4*MassB*Conjg(MassB)/75._dp + 256._dp*g12*g32*MassB*Conjg(MassB)& 
& /9._dp + 128._dp*g12*g32*MassG*Conjg(MassB)/9._dp + 128._dp*g12*g32*MassB*Conjg(MassG)& 
& /9._dp + 256._dp*g12*g32*MassG*Conjg(MassG)/9._dp + 2320._dp*g3**4*MassG*Conjg(MassG)& 
& /3._dp + 32._dp*g1**4*Tr2(1)/15._dp + 80._dp*g3**4*Tr2(3)/3._dp + 16._dp*g12*Tr3(1)  & 
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
&  10._dp*TrYzml2adjYzYzadjYz - TrYzml2CYtYtadjYz + 409._dp*g1**4*MassB*Conjg(MassB)     & 
& /75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)/5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)& 
& /45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)/45._dp + (g12*g22*MassWB*Conjg(MassB))& 
& /5._dp - 4._dp*g12*TrAYzadjYz*Conjg(MassB)/5._dp + 8._dp*g12*MassB*TrYzadjYz*Conjg(MassB)& 
& /5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)/45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)& 
& /45._dp + 32._dp*g22*g32*MassG*Conjg(MassG) + 544._dp*g3**4*MassG*Conjg(MassG)     & 
& /3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG) + (g12*g22*MassB*Conjg(MassWB))    & 
& /5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB) + 2._dp*g12*g22*MassWB*Conjg(MassWB)& 
& /5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB) + 32._dp*g22*g32*MassWB*Conjg(MassWB)  & 
&  + 2._dp*g1**4*Tr2(1)/15._dp + 6._dp*g2**4*Tr2(2) + 32._dp*g3**4*Tr2(3)/3._dp +        & 
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
betamzb22 = 409._dp*g1**4*MassB*Conjg(MassB)/75._dp + 2._dp*g12*g22*MassB*Conjg(MassB)& 
& /5._dp + 32._dp*g12*g32*MassB*Conjg(MassB)/45._dp + 16._dp*g12*g32*MassG*Conjg(MassB)& 
& /45._dp + (g12*g22*MassWB*Conjg(MassB))/5._dp + 16._dp*g12*g32*MassB*Conjg(MassG)& 
& /45._dp + 32._dp*g12*g32*MassG*Conjg(MassG)/45._dp + 32._dp*g22*g32*MassG*Conjg(MassG)& 
&  + 544._dp*g3**4*MassG*Conjg(MassG)/3._dp + 16._dp*g22*g32*MassWB*Conjg(MassG)     & 
&  + (g12*g22*MassB*Conjg(MassWB))/5._dp + 16._dp*g22*g32*MassG*Conjg(MassWB)    & 
&  + 2._dp*g12*g22*MassWB*Conjg(MassWB)/5._dp + 159._dp*g2**4*MassWB*Conjg(MassWB)   & 
&  + 32._dp*g22*g32*MassWB*Conjg(MassWB) + 2._dp*g1**4*Tr2(1)/15._dp +               & 
&  6._dp*g2**4*Tr2(2) + 32._dp*g3**4*Tr2(3)/3._dp - 4._dp*g12*Tr3(1)/3._dp -           & 
&  4._dp*g32*Tr3(3) - 4._dp*(g32*Tr3(3))/sqrt3

 
Dmzb2 = oo16pi2*( betamzb21 + oo16pi2 * betamzb22 ) 

 
Else 
Dmzb2 = oo16pi2* betamzb21 
End If 
 
 
!-------------------- 
! MassB 
!-------------------- 
 
betaMassB1 = 136._dp*g12*MassB/5._dp

 
If (TwoLoopRGE) Then 
betaMassB2 = 6008._dp*g1**4*MassB/75._dp + 348._dp*g12*g22*MassB/5._dp +          & 
&  368._dp*g12*g32*MassB/3._dp + 368._dp*g12*g32*MassG/3._dp + 348._dp*g12*g22*MassWB/5._dp +& 
&  28._dp*g12*TrAYdadjYd/5._dp + 36._dp*g12*TrAYeadjYe/5._dp + 52._dp*g12*TrAYuadjYu/5._dp +& 
&  28._dp*g12*TrAYzadjYz/5._dp + 48._dp*g12*TrCYsAYs/5._dp - 42._dp*g12*MassB*TrCYsYs/5._dp +& 
&  54._dp*g12*TrCYtAYt/5._dp - 9._dp*g12*MassB*TrCYtYt - 28._dp*g12*MassB*TrYdadjYd/5._dp -& 
&  36._dp*g12*MassB*TrYeadjYe/5._dp - 6._dp*g12*MassB*TrYsCYs/5._dp - 9._dp*g12*MassB*TrYtCYt/5._dp -& 
&  52._dp*g12*MassB*TrYuadjYu/5._dp - 28._dp*g12*MassB*TrYzadjYz/5._dp -             & 
&  54._dp*g12*L1*MassB*Conjg(L1)/5._dp + 54._dp*g12*AL1*Conjg(L1)/5._dp -            & 
&  54._dp*g12*L2*MassB*Conjg(L2)/5._dp + 54._dp*g12*AL2*Conjg(L2)/5._dp

 
DMassB = oo16pi2*( betaMassB1 + oo16pi2 * betaMassB2 ) 

 
Else 
DMassB = oo16pi2* betaMassB1 
End If 
 
 
!-------------------- 
! MassWB 
!-------------------- 
 
betaMassWB1 = 16._dp*g22*MassWB

 
If (TwoLoopRGE) Then 
betaMassWB2 = 116._dp*g12*g22*MassB/5._dp + 80._dp*g22*g32*MassG +            & 
&  116._dp*g12*g22*MassWB/5._dp + 376._dp*g2**4*MassWB + 80._dp*g22*g32*MassWB + & 
&  12._dp*g22*TrAYdadjYd + 4._dp*g22*TrAYeadjYe + 12._dp*g22*TrAYuadjYu +          & 
&  12._dp*g22*TrAYzadjYz + 14._dp*g22*TrCYtAYt - 35._dp*g22*MassWB*TrCYtYt/3._dp - & 
&  12._dp*g22*MassWB*TrYdadjYd - 4._dp*g22*MassWB*TrYeadjYe - 7._dp*g22*MassWB*TrYtCYt/3._dp -& 
&  12._dp*g22*MassWB*TrYuadjYu - 12._dp*g22*MassWB*TrYzadjYz - 14._dp*g22*L1*MassWB*Conjg(L1)& 
&  + 14._dp*g22*AL1*Conjg(L1) - 14._dp*g22*L2*MassWB*Conjg(L2) + 14._dp*g22*AL2*Conjg(L2)

 
DMassWB = oo16pi2*( betaMassWB1 + oo16pi2 * betaMassWB2 ) 

 
Else 
DMassWB = oo16pi2* betaMassWB1 
End If 
 
 
!-------------------- 
! MassG 
!-------------------- 
 
betaMassG1 = 8._dp*g32*MassG

 
If (TwoLoopRGE) Then 
betaMassG2 = 46._dp*g12*g32*MassB/3._dp + 46._dp*g12*g32*MassG/3._dp +        & 
&  30._dp*g22*g32*MassG + 1600._dp*g3**4*MassG/3._dp + 30._dp*g22*g32*MassWB +   & 
&  8._dp*g32*TrAYdadjYd + 8._dp*g32*TrAYuadjYu + 8._dp*g32*TrAYzadjYz +            & 
&  18._dp*g32*TrCYsAYs - 63._dp*g32*MassG*TrCYsYs/4._dp - 8._dp*g32*MassG*TrYdadjYd -& 
&  9._dp*g32*MassG*TrYsCYs/4._dp - 8._dp*g32*MassG*TrYuadjYu - 8._dp*g32*MassG*TrYzadjYz

 
DMassG = oo16pi2*( betaMassG1 + oo16pi2 * betaMassG2 ) 

 
Else 
DMassG = oo16pi2* betaMassG1 
End If 
 

 
Call ParametersToG353(Dg1,Dg2,Dg3,DYu,DYd,DYe,DYt,DYs,DYz,DL1,DL2,DMTM,               & 
& Dmue,DMZM,DMSM,DAYu,DAYd,DAYe,DAYt,DAYs,DAYz,DAL1,DAL2,DAmue,DAMTM,DAMZM,              & 
& DAMSM,Dmq2,Dml2,DmHd2,DmHu2,Dmd2,Dmu2,Dme2,Dmt2,Dmtb2,Dms2,Dmsb2,Dmz2,Dmzb2,           & 
& DMassB,DMassWB,DMassG,f)





Iname = Iname - 1 
 
End Subroutine rge353

#endif SARAH


 Subroutine Set_Decoupling_Heavy_States(set)
  Implicit none
  Logical, Intent(in) :: set

   decoupling_heavy_states = set

 End Subroutine Set_Decoupling_Heavy_States


End Module RGEs

