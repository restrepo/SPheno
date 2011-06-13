Module StandardModel

! load modules
Use Control
Use Mathematics
Use RGEs
! load modules


! global variables
 ! Z-boson
 Real(dp), Save :: mZ,mZ2,gamZ,gamZ2,gmZ,gmZ2,BrZqq(5),BrZll(3),BrZinv 
 ! W-boson
 Real(dp), Save :: mW,mW2,gamW,gamW2,gmW,gmW2,BrWqq(2),BrWln(3)
 ! fermion masses
 Real(dp), Save :: mf_l(3), mf_l2(3), mf_u(3), mf_u2(3), mf_d(3) &
   & , mf_d2(3), mf_nu(3), mf_l_mZ(3), mf_d_mZ(3), mf_u_mZ(3)
 ! scale at which the light quark masses for u, d, c, s are defined
 Real(dp), Save :: Q_light_quarks
 ! fermion decay widths
 Real(dp), Save :: g_T = 1.5_dp
 ! couplings
 Real(dp), Save :: Alpha, Alpha_mZ, AlphaS_mZ, G_F, Alpha_mZ_MS, AlphaS_mB
 Real(dp), Save :: Delta_Alpha_Lepton, Delta_Alpha_Hadron
! pseudo observables
 Real(dp), Save :: Rho_parameter
 ! partial widht of mu and tau
 Real(dp), Save :: GammaMu, GammaTau

 Real(dp), Save :: KFactorLee  ! for ISR corrections in e+ e- 

 Complex(dp), Save :: CKM(3,3)
 Real(dp), Save :: lam_wolf=0.2265_dp, A_wolf=0.807_dp, rho_wolf=0.141_dp &
                 & , eta_wolf=0.343_dp

 Complex(dp), Save :: Unu(3,3)
 Real(dp), Save :: theta_12, theta_13, theta_23, delta_nu, alpha_nu1, alpha_nu2
 ! h-bar
 Real(dp), Save :: hbar = 6.58211889e-25_dp
! global variables


Contains

 Subroutine CalculateRunningMasses(mf_l_in, mf_d_in, mf_u_in, Qlow, alpha &
     &  , alphas, Qhigh, mf_l_out, mf_d_out, mf_u_out, kont)
 !-----------------------------------------------------------------------
 ! calculates running masses except the top in the MSbar scheme
 ! at the scale Qout. The formulas for the decoupling procedure can
 ! be found in K.G.Chetyrkin et al., hep-ph/0004189
 ! input:
 !  mf_l_in ......... lepton onshell masses
 !  mf_d_in ......... d-quark masses, it is assumed the d and s are given
 !                    at the scale Qlow and that mb(mb)
 !  mf_u_in ......... u-quark masses, it is assumed the u and c are given
 !                    at the scale Qlow 
 !  Qlow ............ low energy scale, must be smaller or equal mb(mb)
 !  Qhigh ........... high energy scale, must be larger than Qlow and mb(mb)
 !  alpha ........... alpha_magnetic at Qhigh
 !  alphas .......... alpha_strong at Qhigh
 ! output:
 !  mf_l_out ........ mf_l(Qhigh)
 !  mf_d_out ........ mf_d(Qhigh)
 !  mf_u_out ........ mf_u(Qhigh) [except top quark ]
 !  kont ............ contains the error, =0 if everything is fine
 ! written by Werner Porod, 28.8.99
 !-----------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) ::  mf_l_in(3), mf_d_in(3), mf_u_in(3), Qlow, alpha &
     &  , alphas, Qhigh
  Real(dp), Intent(out) :: mf_l_out(3), mf_d_out(3), mf_u_out(3)
  Integer, Intent(inout) :: kont

  Real(dp) :: g9(9), g10(10), as_nf, as_nf_minus_1, as_nf_d_pi, aem , tz, dt
  Real(dp), Parameter :: zeta3 = 1.202056903159594285399738161511449990765_dp &
     & , zeta4 = 1.082323233711138191516003696541167902775_dp                 &
     & , B4 = -1.762800087073770864061897634679818807215_dp                   &
     & , c_as(2) = (/ 11._dp / 72._dp                                         &
     &             , 58067._dp / 13824._dp - 82043._dp * zeta3 / 27648._dp /) &
     & , c_m(2) = (/ 89._dp / 432._dp                                         &
     &            , 713._dp/486._dp - B4 / 36._dp - 221._dp * zeta3 / 288._dp &
     &          + 1.25_dp * zeta4  /)
     
  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateRunningMasses"

  kont = 0 
  !-------------------------------------------------------
  ! check if conditions on Qlow and Qhigh are fullfilled
  !-------------------------------------------------------
  If (Qlow.Gt.mf_d_in(3)) Then
   Iname = Iname - 1
   kont = -101
   Call AddError(101)
   Return
  End If
  If ((Qlow.Gt.Qhigh).Or.(Qhigh.Lt.mf_d_in(3))) Then
   Iname = Iname - 1
   kont = - 102
   Call AddError(102)
   Return
  End If

  !-------------------------------------------------------------------------
  ! calculate first coupling at low scales, putting abritary fermion masses
  ! because they are not needed at this stage
  !-------------------------------------------------------------------------
  g10(1) = Sqrt(4._dp * Pi * alphas)
  g10(2) = Sqrt(4._dp * Pi * alpha)
  g10(3:10)= 2._dp

  If (Qhigh.Gt.mf_d_in(3)) Then ! allow that Qhigh==mb(mb)
   tz = Log(mf_d_in(3)/Qhigh) ! at mb(mb) we have to decouple the b-quark
   dt = tz / 50._dp
   Call odeint(g10, 10, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)
  End If
  g9 = g10(1:9)
  !-------------------------------------------
  ! at mb(mb) we have to decouple the b-quark
  !--------------------------------------------
  If (Qlow.Lt.mf_d_in(3)) Then
   as_nf = g9(1)**2 / (4._dp * Pi)
   as_nf_d_pi = as_nf / Pi
   as_nf_minus_1 = as_nf *(1._dp+ as_nf_d_pi**2 *(c_as(1) +c_as(2)*as_nf_d_pi))
   g9(1) = Sqrt(4._dp * Pi * as_nf_minus_1)
   tz = Log(Qlow/mf_d_in(3)) 
   dt = tz / 50._dp
   Call odeint(g9, 9, 0._dp, tz, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)
  End If
  !--------------------------------------------------
  ! lepton masses at Qlow, note that aem is alpha/pi
  !--------------------------------------------------
  aem = g9(2)**2 / (4._dp * Pi**2)
  g9(3:5) = mf_l_in * (1._dp - aem)
  !-----------
  ! m_u, m_c
  !-----------
  g9(6:7) = mf_u_in(1:2)
  !-----------
  ! m_d, m_s
  !-----------
  g9(8:9) = mf_d_in(1:2)

  !-------------------------
  ! running back to mb(mb)
  !-------------------------
  If (Qlow.Lt.mf_d_in(3)) Then
   tz = Log(Qlow/mf_d_in(3))
   dt = - tz / 50._dp
   Call odeint(g9, 9, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)
   !-------------
   ! thresholds
   !-------------
   AlphaS_mB = as_nf
   g9(1) = Sqrt(4._dp * Pi * as_nf)
   g9(6:9) = g9(6:9) / (1._dp + as_nf_d_Pi**2 * (c_m(1) + c_m(2)*as_nf_d_Pi))  
  End If
  g10(1:9) = g9
  g10(10) = mf_d_in(3)
  !-------------------------
  ! running back to Qhigh
  !-------------------------
  tz = Log(mf_d_in(3)/Qhigh)
  dt = - tz / 50._dp
  Call odeint(g10, 10, tz, 0._dp, 1.e-7_dp, dt, 0._dp, RGE10_SM, kont)

  mf_l_out = g10(3:5)
  mf_u_out(1:2) = g10(6:7)
  mf_u_out(3) = mf_u_in(3)
  mf_d_out = g10(8:10)

  Iname = Iname - 1

 End Subroutine CalculateRunningMasses

 Subroutine FermionMass(Yuk,vev,mF,U_L,U_R,kont)
 !-----------------------------------------------------------------
 ! calculates fermion masses for a given 3*3 yukawa matrix and a vev
 ! input
 !  Yuk ....... 3*3 Yukawa matrix
 !  vev ....... veve
 ! output
 !  mF(i) ..... fermion masses
 !  U_L(i,j) .. left mixing matrix
 !  U_R(i,j) .. right mixing matrix
 ! written by Werner Porod, 30.07.01
 !----------------------------------------------------------------------
 Implicit None

  Complex(Dp), Intent(in) :: Yuk(3,3)
  Real(Dp), Intent(in) :: vev
  Complex(Dp), Intent(out) :: U_L(3,3), U_R(3,3)
  Real(Dp), Intent(out) :: mF(3)
  Integer, Intent(inout) :: kont

  Complex(Dp) :: mat3(3,3)=0, mat32(3,3)=0, U3(3,3), V3(3,3), phaseM
  Real(dp) :: v1(3,3), u1(3,3), mC2(3), test(2)
  Integer :: ierr, i1

  Iname = Iname + 1
  NameOfUnit(Iname) = "FermionMass"

  mF = 0
  U_L = 0
  U_R = 0

  mat3 = oosqrt2 * vev * Yuk
  mat32 = mat3

  mat32 = Matmul( Transpose( Conjg( mat3 ) ), mat3 )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, v1, ierr, test)
   v3 = v1
  Else
   Call EigenSystem(mat32, mC2, v3, ierr, test)
  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( mat3, Transpose( Conjg( mat3 ) ) )
  If ( Maxval( Abs( Aimag(mat32) ) ).Eq.0._dp) Then ! reel matrix
   Call EigenSystem(Real(mat32,dp), mC2, u1, ierr, test)
   u3 = u1
  Else
   Call EigenSystem(mat32, mC2, u3, ierr, test)
  End If
  u3 = Conjg(u3)

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  mat32 = Matmul( Matmul( Conjg(u3), mat3), Transpose( Conjg(v3) ) )
  Do i1=1,3
   phaseM =   mat32(i1,i1)   / Abs( mat32(i1,i1) )
   v3(i1,:) = phaseM * v3(i1,:)
  End Do

  Do i1=1,3
   phaseM = v3(i1,i1)   / Abs( v3(i1,i1) )
   v3(i1,:) = Conjg(phaseM) * v3(i1,:)
   u3(i1,:) = phaseM * u3(i1,:)
  End Do

  If (ierr.Ne.0) Then
   Write(ErrCan,*) 'Warning in subroutine FermionMass, ierr = ',ierr
   Write(ErrCan,*) 'Yuk ',Yuk
   Write(ErrCan,*) 'vev ',vev
   Write(ErrCan,*) ' '
   kont = ierr
   mF(1) = Abs(yuk(1,1)*vev)
   mF(2) = Abs(yuk(2,2)*vev)
   mF(3) = Abs(yuk(3,3)*vev)
   U_L = id3C
   U_R = U_L
   Iname = Iname - 1
   Return
  Endif

  Do i1=1,2
   If ((mC2(i1).Lt.0._dp).And.(Abs(mc2(i1)/mc2(3)).Lt.Epsilon(1._dp))) &
          &  mc2(i1) = 0._dp
  End Do

  mF = Sqrt(MC2) 
  U_L = u3
  U_R = v3

  Iname = Iname - 1

 End Subroutine FermionMass


 Subroutine InitializeStandardModel
 !-----------------------------------------------------------------
 ! reads in Standard Model Paramters from the file StandardModel.in
 ! In case that the file does not exist, default values are used
 ! defined after label 200
 ! written by Werner Porod, 07.07.02
 !------------------------------------------------------------------
 Implicit None

  Integer          :: i1, kont
  Real(dp) :: oo_alpha, s12, s23, s13, phase, c13, c23, c12 &
                &, TimeMu, TimeTau 

  Iname = Iname + 1
  NameOfUnit(Iname) = "InitializeStandardModel"

  !------------------------------------------------------------------------
  ! Contributions to alpha(m_Z), based on F. Jegerlehner, hep-ph/0310234
  ! and Fanchiotti, Kniehl, Sirlin PRD 48 (1993) 307
  !------------------------------------------------------------------------
  Delta_Alpha_Lepton = 0.04020_dp
  Delta_Alpha_Hadron = 0.027651_dp

  Open(99,file='StandardModel.in',status='old',err=200)

 !---------
 ! Z-boson  
 !---------
  Read (99,800) mZ      ! mass
  Read (99,800) gamZ    ! width
  Read (99,*) BrZqq     ! branching ratios in q \bar{q} 
  Read (99,*) BrZll     ! branching ratios in leptons
  Read (99,800) BrZinv  ! invisible branching ratio

  mZ2 = mZ**2
  gamZ2 = gamZ**2
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2
 !---------
 ! W-boson
 !---------
  Read (99,800) mW      ! mass
  Read (99,800) gamW    ! width
  Read (99,*) BrWqq     ! branching ratios in q \bar{q}
  Read (99,*) BrWln     ! branching ratios in leptons

  mW2 = mW**2
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2

 !-----------------------------
 ! lepton masses: e, muon, tau
 !-----------------------------
  Read (99,800) mf_l(1)
  Read (99,800) mf_l(2)
  Read (99,800) mf_l(3)

 !---------------------------------------------------------
 ! scale where masses of light quarks are defined [in GeV]
 !---------------------------------------------------------
  Read (99,*) Q_light_quarks
 !--------------------------
 ! up-quark masses: u, c, t
 !--------------------------
  Read (99,800) mf_u(1)
  Read (99,800) mf_u(2)
  Read (99,800) mf_u(3)

 !----------------------------
 ! down-quark masses: d, s, b
 !----------------------------
  Read (99,800) mf_d(1)
  Read (99,800) mf_d(2)
  Read (99,800) mf_d(3)

  Do i1=1,3
   mf_l2(i1) = mf_l(i1)**2
   mf_u2(i1) = mf_u(i1)**2
   mf_d2(i1) = mf_d(i1)**2
  Enddo
 !--------------------------------------------------------------------
 ! couplings: Alpha(Q=0), Alpha(mZ), Alpha_S(mZ), Fermi constant G_F
 !--------------------------------------------------------------------
  Read (99,800) oo_alpha
  Alpha = 1._dp / oo_alpha

  Read (99,800) oo_alpha
  Alpha_mZ = 1._dp / oo_alpha

  Read (99,800) AlphaS_mZ
  Read (99,800) G_F

 !-----------------------------------------
 ! for ISR correction in e+ e- annihilation
 !-----------------------------------------
  KFactorLee = 1._dp + (Pi/3._dp - 1/(2._dp* Pi) ) * Alpha

 !------------
 ! CKM matrix
 !------------
  Read (99,800) s12
  Read (99,800) s23
  Read (99,800) s13
  Read (99,800) phase

  c12 = Sqrt(1._dp-s12*s12)
  c23 = Sqrt(1._dp-s23*s23)
  c13 = Sqrt(1._dp-s13*s13)

  CKM(1,1) = c12 * c13
  CKM(1,2) = s12 * c13
  CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
  CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,3) = s23 * c13
  CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,3) = c23 * c13

  !-----------------------------
  ! width of mu- and tau-lepton
  !-----------------------------
  Read (99,800) TimeMu
  Read (99,800) TimeTau
  
  GammaMu = G_F**2 * mf_l(2)**5 * (1._dp-8._dp * mf_l(1)**2 / mf_l(2)**2 ) &
        & * (1._dp + 0.5_dp * Alpha * (6.25-Pi2)/Pi) / (192._dp*pi*pi2)
  GammaTau = hbar / TimeTau

  Close(99)

  Call CalculateRunningMasses(mf_l, mf_d, mf_u, Q_light_quarks, alpha_mZ &
     &  , alphas_mZ, mZ, mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)

  Iname = Iname - 1
  Return
!  200 Write(*,*) "File StandardModel.in does not exist, using default values."
!  Write(ErrCan,*) "File StandardModel.in does not exist, using default values."

 !---------
 ! Z-boson  
 !---------
200  mZ = 91.187_dp            ! mass
  gamZ = 2.49_dp            ! width
  BrZqq = 0.14_dp           ! branching ratios in q \bar{q} 
  BrZll(1:2) = 0.035_dp     ! branching ratios in leptons
  BrZll(3) = 0.03_dp 
  BrZinv = 0.2_dp ! invisible branching ratio

  mZ2 = mZ**2
  gamZ2 = gamZ**2
  gmZ = gamZ * mZ
  gmZ2 = gmZ**2
 !---------
 ! W-boson
 !---------
  mW = 80.40_dp      ! mass
  gamW = 2.06_dp     ! width
  BrWqq = 0.35_dp    ! branching ratios in q \bar{q}
  BrWln = 0.1_dp     ! branching ratios in leptons

  mW2 = mW**2
  gamW2 = gamW**2
  gmW = gamW * mW
  gmW2 = gmW**2

 !-----------------------------
 ! lepton masses: e, muon, tau
 !-----------------------------
  mf_l(1) = 0.51099891e-3_dp
  mf_l(2) = 0.105658_dp
  mf_l(3) = 1.7768_dp

 !---------------------------------------------------------
 ! scale where masses of light quarks are defined [in GeV]
 !---------------------------------------------------------
  Q_light_quarks = 2._dp
 !--------------------------
 ! up-quark masses: u, c, t
 !--------------------------
  mf_u(1) = 0.003_dp 
  mf_u(2) = 1.2_dp
  mf_u(3) = 171.3_dp

 !----------------------------
 ! down-quark masses: d, s, b
 !----------------------------
  mf_d(1) = 0.006_dp
  mf_d(2) = 0.105_dp
  mf_d(3) = 4.2_dp

  Do i1=1,3
   mf_l2(i1) = mf_l(i1)**2
   mf_u2(i1) = mf_u(i1)**2
   mf_d2(i1) = mf_d(i1)**2
  Enddo

 !--------------------------------------------------------------------
 ! couplings: Alpha(Q=0), Alpha(mZ), Alpha_S(mZ), Fermi constant G_F
 !--------------------------------------------------------------------
  Alpha = 1._dp / 137.0359895_dp
  Alpha_mZ = 1._dp / 127.9_dp

  AlphaS_mZ = 0.119_dp
  G_F = 1.16639e-5_dp

 !-----------------------------------------
 ! for ISR correction in e+ e- annihilation
 !-----------------------------------------
  KFactorLee = 1._dp + (Pi/3._dp - 1/(2._dp* Pi) ) * Alpha

 !------------
 ! CKM matrix
 !------------
   s12 = lam_wolf
   s23 = s12**2 * A_wolf
   s13 = s23 * lam_wolf * Sqrt(eta_wolf**2+rho_wolf**2) 
   phase = Atan(eta_wolf/rho_wolf)

  c12 = Sqrt(1._dp-s12*s12)
  c23 = Sqrt(1._dp-s23*s23)
  c13 = Sqrt(1._dp-s13*s13)

  CKM(1,1) = c12 * c13
  CKM(1,2) = s12 * c13
  CKM(1,3) = s13 * Exp( (0._dp,-1._dp) * phase )
  CKM(2,1) = -s12*c23 -c12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,2) = c12*c23 -s12*s23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(2,3) = s23 * c13
  CKM(3,1) = s12*s23 -c12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,2) = -c12*s23 - s12*c23*s13 * Exp( (0._dp,1._dp) * phase )
  CKM(3,3) = c23 * c13

  !-----------------------------
  ! width of mu- and tau-lepton
  !-----------------------------
  TimeMu = 2.19709e-6_dp
  TimeTau = 2.9e-13_dp
  
  GammaMu = G_F**2 * mf_l(2)**5 * (1._dp-8._dp * mf_l(1)**2 / mf_l(2)**2 ) &
        & * (1._dp + 0.5_dp * Alpha * (6.25-Pi2)/Pi) / (192._dp*pi*pi2)
  GammaTau = hbar / TimeTau

  Call CalculateRunningMasses(mf_l, mf_d, mf_u, Q_light_quarks, alpha_mZ &
     &  , alphas_mZ, mZ, mf_l_mZ, mf_d_mZ, mf_u_mZ, kont)

  800  Format(f16.7)

  Iname = Iname - 1

 End Subroutine InitializeStandardModel


 Subroutine NeutrinoMasses(MnuL5, mN, N, kont)
 !----------------------------------------------------------------------
 ! calculates neutrino masses + mixing matrix N from the dim 5 operator
 ! input:
 !  MnuL5 .......... dim 5 operator
 ! output 
 !  mN(i) .......... neutrino mass_i
 !  N(i,j) ......... neutrino mixing matrix
 ! written by Werner Porod, 09.01.2006
 ! Note the factor of 1/2 appears because in the RGE running the Higgs field
 ! is still the complex field
 !----------------------------------------------------------------------
 Implicit None

  Integer, Intent(inout) :: kont
  Complex(Dp), Intent(in) :: MnuL5(3,3)
  Real(Dp), Intent(out) :: mN(3)
  Complex(Dp), Intent(out) :: N(3,3)

  Integer :: i1,i2,ierr
  Complex(Dp) :: mat32(3,3), E3(3), phaseM, mat3(3,3)
  Real(Dp) :: g1, gp1, N3a(3,3), eig(3), test(2)


  Iname = Iname + 1
  NameOfUnit(Iname) = 'NeutrinoMasses'

  mat3 = 0.5_dp * MnuL5

  If (Maxval(Abs(Aimag(Mat3))).Lt.  &                           ! matrix is reel
     & (Maxval(Abs(Real(Mat3,dp)))*1.e-3_dp * Epsilon(1._dp))) Then

   Call EigenSystem(Real(Mat3,dp), Eig, N3a, ierr, test)

   Do i1=1,3
    If (Eig(i1).Lt.0._dp) Then
     mN(i1) = - Eig(i1)
     N(i1,:) = (0._dp,1._dp) * N3a(i1,:)
    Else
     mN(i1) = Eig(i1)
     N(i1,:) =N3a(i1,:)
    End If
   End Do

   Do i1=1,2
    Do i2=i1+1,3
     If (mN(i1).Gt.mN(i2)) Then
      Eig(1) = mN(i1)
      mN(i1) = mN(i2)
      mN(i2) = Eig(1)
      E3 = N(i1,:)
      N(i1,:) = N(i2,:)
      N(i2,:) = E3
     End If
    End Do
   End Do

  Else

   mat32 = Matmul( Transpose(Conjg( Mat3 ) ), Mat3 )
   Call EigenSystem(mat32, Eig, N, ierr, test)
   mat32 = Matmul(Conjg(N), Matmul( Mat3, Transpose( Conjg( N ) ) ) )
   Do i1=1,3
    phaseM =   Sqrt( mat32(i1,i1)   / Abs( mat32(i1,i1) ) )
    N(i1,:) = phaseM * N(i1,:)
   End Do
   mN = Sqrt( Eig )

  End If

  If ((ierr.Eq.-14).Or.(ierr.Eq.-16)) Then
    Write(ErrCan,*) "Possible numerical problem in "//NameOfUnit(Iname)
    Write(ErrCan,*) "test =",test
    Write(ErrCan,*) " "
    If (ErrorLevel.Eq.2) Call TerminateProgram
    ierr = 0
  End If
  
  If (ierr.Ne.0) Then
   Write (Errcan,*) 'Warning in subroutine NeutrinoMasses, ierr =',ierr
   Write (Errcan,*) 'MnuL5:',Cmplx(MnuL5(1,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(2,:))
   Write (Errcan,*) '      ',Cmplx(MnuL5(3,:))
   kont = ierr
   Iname = Iname - 1
   Return
  Endif
  !-------------------------------------------
  ! standard convention in neutrino physics
  !-------------------------------------------
  N = Transpose( Conjg( N ) )

  Iname = Iname - 1

 End Subroutine NeutrinoMasses

End Module StandardModel
