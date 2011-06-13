Module Model_Data
!-------------------------------------------------------------
! This module contains global definition of MSSM parameters,
! SUSY and Higgs masses, mixing matrices, decay widths and
! production cross sections.
! The SUSY parameters (denoted generically xxx) are given at four scales:
! mZ ............... xxx_mZ
! Q_EWSB ........... xxx
! M_GUT, M_GMSB .... xxx_0
! M_nu_R ........... xxx_mR
! Some of them may not be given at all scales
!-------------------------------------------------------------
Use Control
Use StandardModel
!------------------------------
! parameters of the Lagrangian
!-----------------------------
 !------------------
 ! gauge couplings
 !------------------
 Real(dp), Dimension(3) :: gauge, gauge_mZ, gauge_0, gauge_mR, gauge_mH3
 !----------------------------
 ! Yukawa couplings
 !----------------------------
 Complex(dp), Dimension(3,3) :: Y_l, Y_d, Y_u, Y_l_mZ, Y_d_mZ, Y_u_mZ &
     & , Y_l_0, Y_nu_0, Y_d_0, Y_u_0, Y_nu, Y_T_0, Y_T_MH3, Y_l_mH3   &
     & , Y_d_mH3, Y_u_mH3, Y_T, Y_S_MH3, Y_Z_MH3
 ! there are 3 R-neutrinos, first index is nu_R index
 Complex(dp), Dimension(3,3,3) :: Y_l_mR, Y_nu_mR, Y_d_mR, Y_u_mR
 Complex(dp) :: lam_0, lamp_0, lam12_0(2), lam12_mH3(2), lam12(2)
 !----------------------------
 ! bilinear Higgs parameters
 !----------------------------
 Complex(dp) :: phase_mu, mu, B, mu_mZ, B_mZ, mu_mR, B_mR, mu_0, B_0, B_MH3 &
     & , mu_MH3, M_RS, M_phi, BM_rs, BM_phi
 !----------------------------
 ! bilinear R-parity parameters
 !----------------------------
 Complex(dp), Dimension(3) :: eps, Beps
 Real(dp), Dimension(3) :: vevL, lam_ex
 !----------------------------
 ! gaugino mass parameters
 !----------------------------
 Complex(dp), Dimension(3) :: Mi, Mi_mZ, Mi_mR, Mi_0, Mi_mH3
 !----------------------------
 ! trilinear couplings
 !----------------------------
 Complex(dp), Dimension(3,3) :: A_l, A_d, A_u, A_l_mZ, A_d_mZ, A_u_mZ, A_l_0   &
     & , A_nu_0, A_d_0, A_u_0, A_nu, A_T_0, A_l_MH3, A_T_MH3, A_d_MH3, A_u_MH3 &
     & , A_T, A_S_MH3, A_Z_MH3
 Complex(dp), Dimension(3,3,3) :: A_l_mR, A_nu_mR, A_d_mR, A_u_mR
 Complex(dp) :: Alam_0, Alamp_0, Alam12_0(2), Alam12_mH3(2), Alam12(2)
 !------------------------------------------------------
 ! trilinear couplings divided by the Yukawa couplings
 !------------------------------------------------------
 Complex(dp), Dimension(3,3) :: AoY_l_0, AoY_nu_0, AoY_d_0, AoY_u_0, AoY_q_0   &
       & , AoY_l, AoY_nu, AoY_d, AoY_u, AoY_l_mZ, AoY_nu_mZ, AoY_d_mZ          &
       & , AoY_u_mZ, AoT_0, AoT_MH3
 Complex(dp) :: Aolam_0, Apolamp_0, Aolam12_0(2), Aolam12_mH3(2)
 !------------------------------------
 ! trilinear couplings, NMSSM
 ! Ao_h0 = A_h0/h0, Ao_lam=A_lam/lam
 ! + spontaneous R-parity violation
 !------------------------------------
 Complex(dp) :: h0, lam, A_h0, A_lam, Ao_h0, Ao_lam, h_pns, A_pns, Ao_hpns
 !-------------------------------------
 ! singlet vev, NMSSM + spontaneous R-parity violation
 !-------------------------------------
 Real(dp) :: vP, vR, vS
 !----------------------------
 ! trilinear R-parity parameters
 !----------------------------
 Complex(dp), Dimension(3,3,3) :: Rp_lam, Rp_lamp, Rp_Alam, Rp_Alamp &
     & , Rp_AoYlam, Rp_AoYlamp
 !----------------------------------
 ! sfermion mass parameters squared
 !----------------------------------
 Complex(dp), Dimension(3,3) :: M2_E, M2_L, M2_D, M2_U, M2_Q, M2_E_mZ, M2_L_mZ &
     & , M2_D_mZ, M2_U_mZ, M2_Q_mZ, M2_E_0, M2_L_0, M2_R_0, M2_D_0, M2_U_0     &
     & , M2_Q_0, M2_R, M2_E_MH3, M2_L_MH3, M2_D_MH3, M2_U_MH3, M2_Q_MH3
 Complex(dp), Dimension(3,3,3) :: M2_E_mR, M2_R_mR, M2_L_mR, M2_D_mR, M2_U_mR &
     & , M2_Q_mR
 !-----------------------------------------------------------
 ! sfermion mass parameters in the super-CKM basis
 !-----------------------------------------------------------
 Complex(dp), Dimension(3,3) :: M2Q_sckm , M2D_sckm , M2U_sckm, Au_sckm, Ad_sckm
 !----------------------------------------------------------------
 ! Higgs mass parameters, tan(beta) and vacuum expectation values
 !----------------------------------------------------------------
 Real(dp) :: M2_H(2), tanb, vevSM(2), M2_H_mZ(2), tanb_mZ, vevSM_mZ(2)      &
     & , M2_H_mR(2), M2_H_0(2), M2_S_0, M2_T_0(2), M2_T_MH3(2), M2_H_MH3(2) &
     & , M_H3(2), M2_T(2), M2_S_MH3(2), M2_Z_MH3(2), MT15_mH3, MZ15_mH3     &
     & , MS15_mH3, M2_P, M2_S
 Logical :: Fifteen_plet = .False.
 !----------------------------
 ! neutrino dim. 5 operator
 !----------------------------
 Complex(dp), Dimension(3,3) :: MnuL5
 !----------------------------
 ! mass of L- and R-neutrinos 
 !----------------------------
 Real(dp), Dimension(3) :: MnuR
 !-------------------------------------------------------------------
 ! SO(10) models, SO(10) scale and D-term of additional gauge boson
 !-------------------------------------------------------------------
 Real(dp), Save :: M_SO_10=0._dp, D_SO_10=0._dp
 !------------------------------------------------------
 ! Munoz model with 2 NuR, additional parameters
 ! Ao_h02 = A_h02/h02, Ao_lam112=A_lam112/lam112
 ! Ao_lam122=A_lam122/lam122, Ao_lam222=A_lam222/lam2
 !------------------------------------------------------
 Real(dp) :: vP2
 Complex(dp) :: h02, lam2, lam112, lam122, A_h02, A_lam112, A_lam122 &
      & , A_lam222, Ao_h02, Ao_lam112, Ao_lam122, Ao_lam222
!------------------------------
! masses and mixing angles
!------------------------------
! MSSM
!----------------------------- 
 !------------------------------------------------------------------
 ! scalar masses (h,H), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS0(2), mS02(2), RS0(2,2), gT_S0(2), gP_S0(2,200), BR_S0(2,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP0(2), mP02(2), RP0(2,2), gT_P0(2), gP_P0(2,200), BR_P0(2,200)
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+), masses squared, corresponding mixing matrix 
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSpm(2), mSpm2(2)
 Complex(dp) :: RSpm(2,2)
 Real(dp) :: gT_Spm(2), gP_Spm(2,200), BR_Spm(2,200)
 !---------------------------------------------------------------------------
 ! chargino masses, masses squared, corresponding mixing matrices
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mC(2), mC2(2)
 Complex(dp) :: U(2,2), V(2,2)
 Real(dp) :: gT_C(2), gP_C(2,267), BR_C(2,267), gP_C2(2,120), BR_C2(2,120) &
     & , BR_C3(2,300), gP_C3(2,300)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) ::  mN(4), mN2(4)
 Complex(dp) :: N(4,4)
 Real(dp) :: gT_N(4), gP_N4_2(4,200), BR_N4_2(4,200) &
           &        , gP_N4_3(4,400), BR_N4_3(4,400)
 Real(dp) :: gT_N5(5), gP_N5(5,350), gP_N(4,350), BR_N(4,350)  &
   & , BR_N5(5,350), gP_N2(5,200), BR_N2(5,200), gP_N3(5,400), BR_N3(5,400)
 !---------------------------------------------------------------------------
 ! gluino mass, phase of the parameter M_3 (=Mi(3))
 ! total decay width, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mGlu
 Complex(dp) :: PhaseGlu
 Real(dp) :: gT_Glu, BR_Glu(230), gP_Glu(230), BR_Glu2(82), gP_Glu2(82) &
      & , BR_Glu3(151), gP_Glu3(151)
 !---------------------------------------------------------------------------
 ! sneutrino masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSneut(3), mSneut2(3)
 Complex(dp) :: Rsneut(3,3)
 Real(dp) :: gT_Sn(3), gP_Sn(3,30), BR_Sn(3,30)
 !---------------------------------------------------------------------------
 ! charged slepton masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSlepton(6), mSlepton2(6)
 Complex(dp) :: RSlepton(6,6)
 Real(dp) :: gT_Sl(6), gP_Sl(6,45), BR_Sl(6,45)
 !---------------------------------------------------------------------------
 ! d-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSdown(6), mSdown2(6)
 Complex(dp) :: RSdown(6,6)
 Real(dp) :: gT_Sd(6), gP_Sd(6,54), BR_Sd(6,54)
 !---------------------------------------------------------------------------
 ! u-squark masses, masses squared, corresponding mixing matrix
 ! total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSup(6), mSup2(6)
 Complex(dp) :: RSup(6,6)
 Real(dp) :: gT_Su(6), gP_Su(6,66), BR_Su(6,66), gP_Su3(6,66), BR_Su3(6,66)
!------------------------------
! NMSSM
!-----------------------------
 !------------------------------------------------------------------
 ! scalar masses (h,H,Re(snu)), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS03(3), mS032(3), RS03(3,3), gT_S03(3)  &
    & , gP_S03(3,200), BR_S03(3,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A,Im(snu)), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP03(3), mP032(3), RP03(3,3), gT_P03(3)  &
    & , gP_P03(3,200), BR_P03(3,200)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 !---------------------------------------------------------------------------
 Real(dp) ::  mN5(5), mN52(5)
 Complex(dp) :: N5(5,5)
  
!------------------------------
! explicit R-parity violation
!----------------------------- 
 !------------------------------------------------------------------
 ! scalar masses (h,H,Re(snu)), masses squared, corresponding mixing matrix, 
 ! total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------
 Real(dp) :: mS05(5), mS052(5), RS05(5,5), gT_S05(5)  &
    & , gP_S05(5,200), BR_S05(5,200)
 !------------------------------------------------------------------------
 ! pseudoscalar masses (G,A,Im(snu)), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !------------------------------------------------------------------------
 Real(dp) :: mP05(5), mP052(5), RP05(5,5), gT_P05(5)  &
    & , gP_P05(5,200), BR_P05(5,200)
 !---------------------------------------------------------------------------
 ! charged scalar masses (G+,H+,sleptons), masses squared, corresponding mixing
 !  matrix, total decay widths, partial decay widths, branching ratios
 !---------------------------------------------------------------------------
 Real(dp) :: mSpm8(8), mSpm82(8)
 Complex(dp) :: RSpm8(8,8)
 Real(dp) :: gT_Spm8(8), gP_Spm8(8,200), BR_Spm8(8,200)
 !---------------------------------------------------------------------------
 ! (lepton,chargino) masses, masses squared, corresponding mixing matrices
 !---------------------------------------------------------------------------
 Real(dp) :: mC5(5), mC52(5)
 Complex(dp) :: U5(5,5), V5(5,5)
 !---------------------------------------------------------------------------
 ! neutralino masses, masses squared, corresponding mixing matrix
 !---------------------------------------------------------------------------
 Real(dp) ::  mN7(7), mN72(7)
 Complex(dp) :: N7(7,7)
 !---------------------------------------------------------------------------
 ! mixing matrices for fermions, mainly needed for later extensions
 !---------------------------------------------------------------------------
 Complex(dp), Dimension(3,3) :: uL_L, uL_R, uD_L, uD_R, uU_L, uU_R
 !-----------------------------------------------------------------------
 ! masses of right handed neutrinos
 !-----------------------------------------------------------------------
 Real(dp) :: m_nu_R(3)

Contains


 Subroutine Set_All_Parameters_0()
 !----------------------------------------------------------------------
 ! this routines sets all model parameters defined above to zero
 !----------------------------------------------------------------------
 Implicit None
  gauge = 0._dp
  gauge_mZ = 0._dp
  gauge_0 = 0._dp
  gauge_mR = 0._dp
  
  Y_l = 0._dp
  Y_d = 0._dp
  Y_u = 0._dp
  Y_l_mZ = 0._dp
  Y_d_mZ = 0._dp
  Y_u_mZ = 0._dp
  Y_l_mR = 0._dp
  Y_nu_mR = 0._dp
  Y_d_mR = 0._dp
  Y_u_mR = 0._dp
  Y_l_0 = 0._dp
  Y_nu_0 = 0._dp
  Y_d_0 = 0._dp
  Y_u_0 = 0._dp
  Y_T_0 = 0._dp
  lam12_0 = 0._dp
  Y_T_MH3 = 0._dp
  Y_Z_MH3 = 0._dp
  Y_S_MH3 = 0._dp
  lam12_MH3 = 0._dp

  phase_mu = 0._dp
  mu = 0._dp
  B = 0._dp
  mu_mZ = 0._dp
  B_mZ = 0._dp
  mu_mR = 0._dp
  B_mR = 0._dp
  mu_0 = 0._dp
  B_0 = 0._dp
  
  eps = 0._dp
  Beps = 0._dp
  vevL = 0._dp
  lam_ex = 0._dp
  
  Rp_lam = 0._dp
  Rp_lamp = 0._dp
  Rp_Alam = 0._dp
  Rp_Alamp = 0._dp
  Rp_AoYlam = 0._dp
  Rp_AoYlamp = 0._dp

  Mi = 0._dp
  Mi_mZ = 0._dp
  Mi_mR = 0._dp
  Mi_0 = 0._dp
  
  A_l = 0._dp
  A_d = 0._dp
  A_u = 0._dp
  A_nu = 0._dp
  A_l_mZ = 0._dp
  A_d_mZ = 0._dp
  A_u_mZ = 0._dp
  A_l_mR = 0._dp
  A_nu_mR = 0._dp
  A_d_mR = 0._dp
  A_u_mR = 0._dp
  A_l_0 = 0._dp
  A_nu_0 = 0._dp
  A_d_0 = 0._dp
  A_u_0 = 0._dp
  AoY_l_0 = 0._dp
  AoY_nu_0 = 0._dp
  AoY_d_0 = 0._dp
  AoY_u_0 = 0._dp
  AoY_q_0 = 0._dp
  AoY_l = 0._dp
  AoY_nu = 0._dp
  AoY_d = 0._dp
  AoY_u = 0._dp
  AoY_l_mZ = 0._dp
  AoY_nu_mZ = 0._dp
  AoY_d_mZ = 0._dp
  AoY_u_mZ = 0._dp
  A_T_0 = 0._dp
  Alam12_0 = 0._dp
  A_T_MH3 = 0._dp
  A_S_MH3 = 0._dp
  A_Z_MH3 = 0._dp
  Alam12_MH3 = 0._dp
  AoT_0 = 0._dp
  Aolam12_0 = 0._dp
  AoT_MH3 = 0._dp
  Aolam12_MH3 = 0._dp
  
  h0 = 0._dp
  lam = 0._dp
  A_h0 = 0._dp
  A_lam = 0._dp
  Ao_h0 = 0._dp
  Ao_lam = 0._dp

  vP = 0._dp
  
  h02 = 0._dp
  lam2 = 0._dp
  lam112 = 0._dp
  lam122 = 0._dp
  A_h02 = 0._dp
  A_lam112 = 0._dp
  A_lam122 = 0._dp
  A_lam222 = 0._dp
  Ao_h02 = 0._dp
  Ao_lam112 = 0._dp
  Ao_lam122 = 0._dp
  Ao_lam222 = 0._dp

  vP2 = 0._dp

  M2_E = 0._dp
  M2_L = 0._dp
  M2_R = 0._dp
  M2_D = 0._dp
  M2_U = 0._dp
  M2_Q = 0._dp
  M2_E_mZ = 0._dp
  M2_L_mZ = 0._dp
  M2_D_mZ = 0._dp
  M2_U_mZ = 0._dp
  M2_Q_mZ = 0._dp
  M2_E_mR = 0._dp
  M2_R_mR = 0._dp
  M2_L_mR = 0._dp
  M2_D_mR = 0._dp
  M2_U_mR = 0._dp
  M2_Q_mR = 0._dp
  M2_E_0 = 0._dp
  M2_L_0 = 0._dp
  M2_R_0 = 0._dp
  M2_D_0 = 0._dp
  M2_U_0 = 0._dp
  M2_Q_0 = 0._dp

  M2Q_sckm = 0._dp
  M2D_sckm = 0._dp
  M2U_sckm = 0._dp
  Au_sckm = 0._dp
  Ad_sckm = 0._dp

  M2_S_MH3 = 0._dp
  M2_Z_MH3 = 0._dp
  MT15_mH3 = 0._dp
  MZ15_mH3 = 0._dp
  MS15_mH3= 0._dp

  M2_H = 0._dp
  tanb = 0._dp
  vevSM = 0._dp
  M2_H_mZ = 0._dp
  tanb_mZ = 0._dp
  vevSM_mZ = 0._dp
  M2_H_mR = 0._dp
  M2_H_0 = 0._dp
  
  MnuL5 = 0._dp
  MnuR = 0._dp

  lam_0 = 0._dp
  Alam_0 = 0._dp
  Aolam_0 = 0._dp
  lamp_0 = 0._dp
  Alamp_0 = 0._dp
  Apolamp_0 = 0._dp

  M_SO_10=0._dp
  D_SO_10=0._dp

  U = 0._dp
  U5 = 0._dp
  V = 0._dp
  V5 = 0._dp
  N = 0._dp
  N5 = 0._dp
  N7 = 0._dp
  RSup = 0._dp
  RSdown = 0._dp
  RSneut = 0._dp
  RSlepton = 0._dp
  
 End Subroutine Set_All_Parameters_0


 Subroutine Switch_from_superCKM(Y_d, Y_u, Ad_in, Au_in, MD_in, MQ_in, MU_in &
                      &, Ad_out, Au_out, MD_out, MQ_out, MU_out, tr        &
                      &, RSd_in, RSu_in, RSd_out, RSu_out, CKM_out, Yd, Yu )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the  super CKM basis to the electroweak basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 implicit none
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Au_in, Ad_in, MD_in &
        & , MQ_in, MU_in
  Complex(dp), Optional, Intent(in), Dimension(6,6) :: RSu_in, RSd_in
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Au_out, Ad_out, MD_out, MQ_out &
        & , MU_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSu_out, RSd_out
  Complex(dp), Optional, Intent(out) :: CKM_out(3,3)
  Real(dp), Optional, Intent(out) :: Yd(3), Yu(3)

  Complex(dp), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R, CKM_Q
  Complex(dp) :: rot(6,6)

  real(dp) :: mf(3)
  integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  if (tr) then
   Call FermionMass(Transpose(Y_u), 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Transpose(Y_d), 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  else
   Call FermionMass(Y_u, 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Y_d, 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  end if
  !---------------------------------------------------------
  ! CKM matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )
  uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  !-------------------------------------------------------------------
  ! also the right quark must be multiplied with the conjugate phase
  ! as otherwise the masses get complex
  !-------------------------------------------------------------------
  uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  uU_R(2,:) = uU_R(2,:) / Conjg(CKM_Q(2,3)) * Abs(CKM_Q(2,3))
  uU_R(3,:) = uU_R(3,:) / Conjg(CKM_Q(3,3)) * Abs(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  If (present(CKM_out)) CKM_out = CKM_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super CKM basis
  !-------------------------------------------------------------------
  if (tr) then
   Au_out = Matmul( Matmul(Transpose(Conjg(uU_R)), Au_in), uU_L)
   Ad_out = Matmul( Matmul(Transpose(Conjg(uD_R)), Ad_in), uD_L)

   MD_out = Matmul( Matmul( Transpose(Conjg(uD_R)), MD_in), uD_R)
   MU_out = Matmul( Matmul( Transpose(Conjg(uU_R)), MU_in), uU_R)
   MQ_out = Matmul( Matmul( Conjg(uD_L), MQ_in), Transpose(uD_L) )

   If (Present(RSu_in).and.Present(RSu_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uU_L
    rot(4:6,4:6) = Conjg(uU_R)
    RSu_out = Matmul(Conjg(rot), RSu_in)
   end if
   If (Present(RSd_in).and.Present(RSd_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uD_L
    rot(4:6,4:6) = Conjg(uD_R)
    RSd_out = Matmul(Conjg(rot), RSd_in)
   end if

  else
   Au_out = Matmul( Matmul(Transpose(uU_L), Au_in), Conjg(uU_R))
   Ad_out = Matmul( Matmul(Transpose(uD_L), Ad_in), Conjg(uD_R))

   MD_out = Matmul( Matmul( Transpose(uD_R), MD_in), Conjg(uD_R))
   MU_out = Matmul( Matmul( Transpose(uU_R), MU_in), Conjg(uU_R))
   MQ_out = Matmul( Matmul( Transpose(Conjg(uD_L)), MQ_in), uD_L )

   If (Present(RSu_in).and.Present(RSu_out)) then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uU_L)
    rot(4:6,4:6) = uU_R
    RSu_out = Matmul(RSu_in, Transpose(rot))
   end if
   If (Present(RSd_in).and.Present(RSd_out)) then
    rot = 0._dp
    rot(1:3,1:3) = Conjg(uD_L)
    rot(4:6,4:6) = uD_R
    RSd_out = Matmul(RSd_in, Transpose(rot))
   end if
  end if


 end Subroutine Switch_from_superCKM

 Subroutine Switch_to_superCKM(Y_d, Y_u, Ad_in, Au_in, MD_in, MQ_in, MU_in &
                      &, Ad_out, Au_out, MD_out, MQ_out, MU_out, tr        &
                      &, RSd_in, RSu_in, RSd_out, RSu_out, CKM_out, Yd, Yu )
 !---------------------------------------------------------------------------
 ! shifts the parameter from the electroweak basis to the super CKM basis
 ! written by werner Porod, 12.03.08
 !---------------------------------------------------------------------------
 implicit none
  Complex(dp), Intent(in), Dimension(3,3) :: Y_d, Y_u, Au_in, Ad_in, MD_in &
        & , MQ_in, MU_in
  Complex(dp), Optional, Intent(in), Dimension(6,6) :: RSu_in, RSd_in
  Logical, Intent(in) :: tr  ! if true, then the matrices are transposed 
                             ! compared to low energy definition
  Complex(dp), Intent(out), Dimension(3,3) :: Au_out, Ad_out, MD_out, MQ_out &
        & , MU_out
  Complex(dp), Optional, Intent(out), Dimension(6,6) :: RSu_out, RSd_out
  Complex(dp), Optional, Intent(out) :: CKM_out(3,3)
  Real(dp), Optional, Intent(out) :: Yd(3), Yu(3)

  Complex(dp), Dimension(3,3) :: uU_L, uU_R, uD_L, uD_R, CKM_Q
  Complex(dp) :: rot(6,6)

  real(dp) :: mf(3)
  integer :: ierr

  !------------------------------------------
  ! diagonalizing d- and u-Yukawa couplings
  ! I am only interested in the mixing matrices
  !------------------------------------------
  if (tr) then
   Call FermionMass(Transpose(Y_u), 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Transpose(Y_d), 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  else
   Call FermionMass(Y_u, 1._dp, mf, uU_L, uU_R, ierr)
   If (Present(Yu)) Yu = sqrt2 * mf
   Call FermionMass(Y_d, 1._dp, mf, uD_L, uD_R, ierr)
   If (Present(Yd)) Yd = sqrt2 * mf
  end if
  !---------------------------------------------------------
  ! CKM matrix at Q, shifting phases according to PDG form
  !---------------------------------------------------------
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )
  uD_L(1,:) = uD_L(1,:) / Conjg(CKM_Q(1,1)) * Abs(CKM_Q(1,1))
  uD_L(2,:) = uD_L(2,:) / Conjg(CKM_Q(1,2)) * Abs(CKM_Q(1,2))
  uU_L(2,:) = uU_L(2,:) / CKM_Q(2,3) * Abs(CKM_Q(2,3))
  uU_L(3,:) = uU_L(3,:) / CKM_Q(3,3) * Abs(CKM_Q(3,3))
  !-------------------------------------------------------------------
  ! also the right quark must be multiplied with the conjugate phase
  ! as otherwise the masses get complex
  !-------------------------------------------------------------------
  uD_R(1,:) = uD_R(1,:) / CKM_Q(1,1) * Abs(CKM_Q(1,1))
  uD_R(2,:) = uD_R(2,:) / CKM_Q(1,2) * Abs(CKM_Q(1,2))
  uU_R(2,:) = uU_R(2,:) / Conjg(CKM_Q(2,3)) * Abs(CKM_Q(2,3))
  uU_R(3,:) = uU_R(3,:) / Conjg(CKM_Q(3,3)) * Abs(CKM_Q(3,3))
  CKM_Q =  Matmul(uU_L, Transpose(Conjg(ud_L)) )

  If (present(CKM_out)) CKM_out = CKM_Q
  !-------------------------------------------------------------------
  ! shifting the parameters to the super CKM basis
  !-------------------------------------------------------------------
  if (tr) then
   Au_out = Matmul( Matmul(uU_R, Au_in), Transpose(Conjg(uU_L)))
   Ad_out = Matmul( Matmul(uD_R, Ad_in), Transpose(Conjg(uD_L)))

   MD_out = Matmul( Matmul( uD_R, MD_in), Transpose(Conjg(uD_R)))
   MU_out = Matmul( Matmul( uU_R, MU_in), Transpose(Conjg(uU_R)))
   MQ_out = Matmul( Matmul( Transpose(uD_L), MQ_in), Conjg(uD_L) )

   If (Present(RSu_in).and.Present(RSu_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uU_L
    rot(4:6,4:6) = Conjg(uU_R)
    RSu_out = Matmul(Transpose(rot), RSu_in)
   end if
   If (Present(RSd_in).and.Present(RSd_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uD_L
    rot(4:6,4:6) = Conjg(uD_R)
    RSd_out = Matmul(Transpose(rot), RSd_in)
   end if

  else
   Au_out = Matmul( Matmul(Conjg(uU_L), Au_in), Transpose(uU_R))
   Ad_out = Matmul( Matmul(Conjg(uD_L), Ad_in), Transpose(uD_R))

   MD_out = Matmul( Matmul( Conjg(uD_R), MD_in), Transpose(uD_R))
   MU_out = Matmul( Matmul( Conjg(uU_R), MU_in), Transpose(uU_R))
   MQ_out = Matmul( Matmul( uD_L, MQ_in), Transpose(Conjg(uD_L)) )

   If (Present(RSu_in).and.Present(RSu_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uU_L
    rot(4:6,4:6) = Conjg(uU_R)
    RSu_out = Matmul(RSu_in, rot)
   end if
   If (Present(RSd_in).and.Present(RSd_out)) then
    rot = 0._dp
    rot(1:3,1:3) = uD_L
    rot(4:6,4:6) = Conjg(uD_R)
    RSd_out = Matmul(RSd_in, rot)
   end if
  end if


 end Subroutine Switch_to_superCKM

End Module Model_Data

