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
Logical, Private, Save :: OnlyDiagonal
! private variables

!------------------------------------------------
! new variables
!------------------------------------------------
 Integer, Private, Save :: loop_ord
 Logical, Private, Save :: Nu_R, Phi, S_f, LH, LLE, UDD, l_mu, l_lin, l_qua

Contains


 Subroutine CalculateAnomalousDimensions(g2, Y_l, Y_d, Y_u, Y_nu, h0, lambda &
              & , h_pnus, LambdaE, LambdaP, LambdaPP  &
              & , Mi, A_l, A_nu, A_d, A_u, Ah0, Alam, A_h &
              & , ME2, MR2, ML2, MD2, MU2, MQ2, MH2, Mphi2, MS2 &
              & , beta, GammaE, GammaR, GammaL, GammaD, GammaU, GammaQ, GammaH &
              & , GammaPhi, GammaS, GammaLH, betaM, GammaAE, GammaAR, GammaAL  &
              & , GammaAD, GammaAU, GammaAQ, GammaAH, GammaAPhi, GammaAS &
              & , DME2, DMR2, DML2, DMD2, DMU2, DMQ2, DMH2, DMphi2, DMS2 )
 !---------------------------------------------------------------------------
 ! Calculates the anomalous dimensions of the chiral superfields in the MSSM
 ! Input: g2 ........... gauge couplings squared
 !        Y_l .......... lepton Yukawa coupling
 !        Y_d .......... d-quark Yukawa coupling
 !        Y_u .......... u-quark Yukawa coupling
 !        Y_nu ......... neutrino Yukawa coupling (optional)
 !        h0 ........... Phi H_1 H_2 coupling (optional)
 !        lambda ....... Phi^3 coupling (optional)
 !        h_pnus ....... Phi Nu^C S coupling (optional)
 !        LambdaE ...... LLE couplig (optional)
 !        LambdaP ...... LQD couplig (optional)
 !        LambdaPP ..... UDD couplig (optional)
 !        Mi ........... gaugino mass parameters
 ! Output: 
 !        beta ......... beta functions of the gauge couplings
 !        GammaE ....... anomalous dimensions of the right sleptons
 !        GammaL ....... anomalous dimensions of the left sleptons
 !        GammaD ....... anomalous dimensions of the right d-squarks
 !        GammaU ....... anomalous dimensions of the right u-squarks
 !        GammaQ ....... anomalous dimensions of the left squarks
 !        GammaH ....... anomalous dimensions of the Higgs fields
 !        GammaR ....... anomalous dimensions of the right neutrinos (optional)
 !        GammaPhi ..... anomalous dimensions of the Phi field (optional)
 !        GammaS ....... anomalous dimensions of the S fields (optional)
 !        betaM ........ beta functions of the gauginos
 !         
 ! written by Werner Porod, 23 March 2004
 ! 26.03.04: - adding right handed neutrinos, 1 + 2 loop
 !           - adding Phi and S fields, 1 loop
 !           - adding trilinear R-parity couplings, 1-loop
 ! 23.02.06: adding beta functions
 ! 24.02.06: - anomalous couplings with dimension [-] and [m] added
 !           - start with dim. [m^2]
 !---------------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: g2(3), Mphi2, MH2(2)
  Complex(dp), Intent(in) :: Mi(3)
  Complex(dp), Intent(in), Dimension(3,3) :: Y_l, Y_d, Y_u, A_l, A_d, A_u &
     & , ME2, ML2, MD2, MU2, MQ2
  Complex(dp), Intent(in), Dimension(3,3), Optional :: Y_nu, A_nu, h_pnus, A_h &
             & , MR2, MS2
  Complex(dp), Intent(in), Dimension(3,3,3), Optional :: LambdaE, LambdaP &
      & , LambdaPP
  Complex(dp), Intent(in), Optional :: h0, lambda, Ah0, Alam

  Real(dp), Intent(out) :: beta(3), GammaH(2), GammaAH(2), DMH2(2), DMphi2
  Complex(dp), Intent(out), Dimension(3) :: GammaLH, BetaM
  Complex(dp), Intent(out), Dimension(3,3) :: GammaE, GammaL, GammaD, GammaU &
        &  , GammaQ
  Complex(dp), Intent(out), Dimension(3,3), Optional :: GammaR, GammaS, GammaAR &
        & , GammaAS, DME2, DMR2, DML2, DMD2, DMU2, DMQ2, DMS2, GammaAE, GammaAL &
        & , GammaAD, GammaAU, GammaAQ
  Complex(dp), Intent(out), Optional :: GammaPhi, GammaAPhi

  Real(dp) :: G_fac(5), TraceY(5), g4(3), sumH(2), TraceY2(6), h02, lambda2 &
        & , g2Mi2(3), Ah02, Alam2, TraceA(5)
  Complex(dp) :: TraceAY(5), g2Mi(3), Ah0ah0, Alam_alam
  Complex(dp), Dimension(3,3) :: Y_la, Y_da, Y_ua, YlaYl, aYlYl, YdaYd, aYdYd  &
       & , YuaYu, aYuYu, sumE, sumL, sumD, sumU, sumQ, YlaYlYlaYl, aYlYlaYlYl  &
       & , YdaYdYdaYd, aYdYdaYdYd, YuaYuYuaYu, aYuYuaYuYu, AlaYl, aYlAl, AdaYd &
       & , aYdAd, AuaYu, aYuAu, AH_AH, A_la, A_da, A_ua, AlaAl, aAlAl, AdaAd   &
       & , aAdAd, AuaAu, aAuAu, A_haA_h, A_hA_ha, A_ha
  Complex(dp), Dimension(3,3) ::Y_nua, YnuaYnu, aYnuYnu, sumR, YnuaYnuYnuaYnu  &
       & , aYnuYnuaYnuYnu, h_pnusA, hah, ahh, AnuaYnu, aYnuAnu, A_nua, AnuaAnu &
       & , aAnuAnu

  Integer :: i1, i2, i3, i4

  Iname = Iname + 1
  NameOfUnit(Iname) = "CalculateAnomalousDimensions"

  !------------------------------------------------------------------------
  ! check which of the additional Yukawas beside MSSM Yukawas are present
  !------------------------------------------------------------------------
 
  GammaE = ZeroC
  GammaL = ZeroC
  GammaD = ZeroC
  GammaU = ZeroC
  GammaQ = ZeroC
  GammaH = 0._dp

  If (Nu_R) GammaR = ZeroC

  If (loop_ord.Gt.2) Then
   Write(ErrCan,*) "Warning from routine "//NameOfUnit(Iname)
   Write(ErrCan,*) "Only 2-loops implemented and not the required ",loop_ord &
                 & ,"loops"
   If (ErrorLevel.Eq.2) Call TerminateProgram
  End If

  OnlyDiagonal = .Not.GenerationMixing

  Call Adjungate(Y_l,Y_la)
  Call Adjungate(Y_d,Y_da)
  Call Adjungate(Y_u,Y_ua)

  YlaYl = MatMul2(Y_l,Y_la,OnlyDiagonal)
  aYlYl = MatMul2(Y_la,Y_l,OnlyDiagonal)
  YdaYd = MatMul2(Y_d,Y_da,OnlyDiagonal)
  aYdYd = MatMul2(Y_da,Y_d,OnlyDiagonal)
  YuaYu = MatMul2(Y_u,Y_ua,OnlyDiagonal)
  aYuYu = MatMul2(Y_ua,Y_u,OnlyDiagonal)

  If (Nu_R) Then
   Call Adjungate(Y_nu,Y_nua)
   YnuaYnu = MatMul2(Y_nu,Y_nua,OnlyDiagonal)
   aYnuYnu = MatMul2(Y_nua,Y_nu,OnlyDiagonal)
  End If

  If (Phi) Then
   h02 = Abs(h0)**2
   lambda2 = Abs(Lambda)**2
  End If

  If (S_f) Then
   Call Adjungate(h_pnus,h_pnusA)
   ahh = MatMul2(h_pnusA,h_pnus, OnlyDiagonal)
   hah = MatMul2(h_pnus,h_pnusA, OnlyDiagonal)
  End If

  If (l_lin) Then
   AlaYl = MatMul2(A_l,Y_la,OnlyDiagonal)
   aYlAl = MatMul2(Y_la,A_l,OnlyDiagonal)
   AdaYd = MatMul2(A_d,Y_da,OnlyDiagonal)
   aYdAd = MatMul2(Y_da,A_d,OnlyDiagonal)
   AuaYu = MatMul2(A_u,Y_ua,OnlyDiagonal)
   aYuAu = MatMul2(Y_ua,A_u,OnlyDiagonal)

   If (Nu_R) Then
    AnuaYnu = MatMul2(A_nu,Y_nua,OnlyDiagonal)
    aYnuAnu = MatMul2(Y_nua,A_nu,OnlyDiagonal)
   End If
   If (Phi) Then
    Ah0ah0 = Ah0 * Conjg(h0)
    Alam_alam = Alam * Conjg(Lambda)
   End If

   If (S_f) Then
    Ah_ah = Matmul(A_h,h_pnusA)
   End If
  End If

  If (l_qua) Then
   Call Adjungate(A_l,A_la)
   Call Adjungate(A_d,A_da)
   Call Adjungate(A_u,A_ua)

   AlaAl = MatMul2(A_l,A_la,OnlyDiagonal)
   aAlAl = MatMul2(A_la,A_l,OnlyDiagonal)
   AdaAd = MatMul2(A_d,A_da,OnlyDiagonal)
   aAdAd = MatMul2(A_da,A_d,OnlyDiagonal)
   AuaAu = MatMul2(A_u,A_ua,OnlyDiagonal)
   aAuAu = MatMul2(A_ua,A_u,OnlyDiagonal)

   If (Nu_R) Then
    Call Adjungate(A_nu,A_nua)
    AnuaAnu = MatMul2(A_nu,A_nua,OnlyDiagonal)
    aAnuAnu = MatMul2(A_nua,A_nu,OnlyDiagonal)
   End If

   If (Phi) Then
    Ah02 = Abs(Ah0)**2
    Alam2 = Abs(Alam)**2
   End If

   If (S_f) Then
    Call Adjungate(A_h,A_ha)
    A_haA_h = MatMul2(A_ha,A_h, OnlyDiagonal)
    A_hA_ha = MatMul2(A_h, A_ha, OnlyDiagonal)
   End If
  End If
  !-------------------------------------
  ! one loop order
  !-------------------------------------
  ! gauge + gauginos
  !------------------
  beta = b_1 * g2
  If (l_lin) betaM = 2._dp * b_1 * g2 * Mi

  !-------------------------------
  ! anomalous dimensions
  !-------------------------------
  GammaE = 2._dp * aYlYl
  GammaL = YlaYl
  GammaD = 2._dp * aYdYd
  GammaU = 2._dp * aYuYu
  GammaQ = YdaYd + YuaYu

  g_fac(1) = 1.2_dp * g2(1)
  g_fac(2) = 0.3_dp * g2(1) + 1.5_dp * g2(2)
  g_fac(3) = 2._dp * g2(1) / 15._dp + 8._dp * g2(3) / 3._dp
  g_fac(4) = 8._dp * g2(1) / 15._dp + 8._dp * g2(3) / 3._dp
  g_fac(5) = 1._dp * g2(1) / 30._dp + 1.5_dp * g2(2) + 8._dp * g2(3) / 3._dp
  Do i1=1,3
   GammaE(i1,i1) = GammaE(i1,i1) - g_fac(1)
   GammaL(i1,i1) = GammaL(i1,i1) - g_fac(2)
   GammaD(i1,i1) = GammaD(i1,i1) - g_fac(3)
   GammaU(i1,i1) = GammaU(i1,i1) - g_fac(4)
   GammaQ(i1,i1) = GammaQ(i1,i1) - g_fac(5)
  End Do
  TraceY(1) = CTrace(YlaYl)
  TraceY(2) = CTrace(YdaYd)
  TraceY(3) = CTrace(YuaYu)
  GammaH(1) = 3._dp * TraceY(2) + TraceY(1) - g_fac(2)   
  GammaH(2) = 3._dp * TraceY(3) - g_fac(2)   

  If (Nu_R) Then
   GammaR = 2._dp * aYnuYnu
   GammaL = GammaL + YnuaYnu
   TraceY(4) = CTrace(YnuaYnu)
   GammaH(2) = GammaH(2) + TraceY(4)
  End If

  If (Phi) Then
   GammaPhi = 2._dp * h02 + 0.5_dp * lambda2
   GammaH = GammaH + h02 
  End If

  If (S_f) Then
   TraceY(5) = CTrace(ahh)
   GammaPhi = GammaPhi + TraceY(5)
   GammaR = GammaR + hah
   GammaS = ahh
  End If

  If (LH) Then
   GammaLH = ZeroC
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      GammaLH(i1) = GammaLH(i1) - 3._dp * Conjg(LambdaP(i1,i2,i3)) * Y_d(i2,i3) &
                  &             -  Conjg(LambdaE(i1,i2,i3)) * Y_l(i2,i3)
      Do i4=1,3
       GammaE(i1,i2) = GammaE(i1,i2) + LambdaE(i3,i4,i1)*Conjg(LambdaE(i3,i4,i2))
       GammaL(i1,i2) = GammaL(i1,i2) &
                   & + LambdaE(i1,i3,i4) * Conjg(LambdaE(i2,i3,i4))         &
                   & + 3._dp * LambdaP(i1,i3,i4) * Conjg(LambdaP(i2,i3,i4))
       GammaD(i1,i2) = GammaD(i1,i2) &
                   & + 2._dp * LambdaP(i3,i4,i2) * Conjg(LambdaP(i3,i4,i1))
       GammaQ(i1,i2) = GammaQ(i1,i2) + LambdaP(i3,i1,i4)*Conjg(LambdaP(i3,i2,i4))
      End Do
     End Do
    End Do
   End Do
  End If

  If (UDD) Then
   Do i1=1,3
    Do i2=1,3
     Do i3=1,3
      Do i4=1,3
       GammaD(i1,i2) = GammaD(i1,i2) &
                   & + 2._dp * LambdaPP(i3,i1,i4) * Conjg(LambdaPP(i3,i2,i4))
       GammaU(i1,i2) = GammaU(i1,i2) &
                   & + LambdaPP(i1,i3,i4) * Conjg(LambdaPP(i2,i3,i4))
      End Do
     End Do
    End Do
   End Do
  End If

  !----------------------------------------------------
  ! anomalous dim., A-parameters
  !----------------------------------------------------
  If (l_lin) Then
   GammaAE = -2._dp * aYlAl
   GammaAL = -AlaYl
   GammaAD = -2._dp * aYdAd
   GammaAU = -2._dp * aYuAu
   GammaAQ = -AdaYd - AuaYu

   g2Mi = g2 * Mi

   g_fac(1) = 1.2_dp * g2Mi(1)
   g_fac(2) = 0.3_dp * g2Mi(1) + 1.5_dp * g2Mi(2)
   g_fac(3) = 2._dp * g2Mi(1) / 15._dp + 8._dp * g2Mi(3) / 3._dp
   g_fac(4) = 8._dp * g2Mi(1) / 15._dp + 8._dp * g2Mi(3) / 3._dp
   g_fac(5) = 1._dp * g2Mi(1) / 30._dp + 1.5_dp * g2Mi(2) &
          & + 8._dp * g2Mi(3) / 3._dp
   Do i1=1,3
    GammaAE(i1,i1) = GammaAE(i1,i1) - g_fac(1)
    GammaAL(i1,i1) = GammaAL(i1,i1) - g_fac(2)
    GammaAD(i1,i1) = GammaAD(i1,i1) - g_fac(3)
    GammaAU(i1,i1) = GammaAU(i1,i1) - g_fac(4)
    GammaAQ(i1,i1) = GammaAQ(i1,i1) - g_fac(5)
   End Do
   TraceAY(1) = CTrace(AlaYl)
   TraceAY(2) = CTrace(AdaYd)
   TraceAY(3) = CTrace(AuaYu)
   GammaAH(1) = -3._dp * TraceAY(2) - TraceAY(1) - g_fac(2)   
   GammaAH(2) = -3._dp * TraceAY(3) - g_fac(2)   

   If (Nu_R) Then
    GammaAR = - 2._dp * aYnuAnu
    GammaAL = GammaAL - AnuaYnu
    TraceAY(4) = CTrace(AnuaYnu)
    GammaAH(2) = GammaAH(2) - TraceAY(4)
   End If

   If (Phi) Then
    GammaAPhi = - 2._dp * Ah0ah0 - 0.5_dp * Alam_alam
    GammaAH = GammaAH - Ah0ah0 
   End If

   If (S_f) Then
    TraceAY(5) = CTrace(Ah_Ah)
    GammaAPhi = GammaAPhi - TraceAY(5)
    GammaAR = GammaAR - Ah_ah
    GammaAS = - Ah_ah
   End If
  End If

  !----------------------------------------------------
  ! anomalous dim., mass squared parameters
  !----------------------------------------------------
  If (l_qua) Then
   DME2 = 2._dp * aAlAl
   DML2 = AlaAl
   DMD2 = 2._dp * aAdAd
   DMU2 = 2._dp * aAuAu
   DMQ2 = AdaAd + AuaAu

   g2Mi2 = g2 * Abs(Mi)**2

   g_fac(1) = 1.2_dp * g2Mi2(1)
   g_fac(2) = 0.3_dp * g2Mi2(1) + 1.5_dp * g2Mi2(2)
   g_fac(3) = 2._dp * g2Mi2(1) / 15._dp + 8._dp * g2Mi2(3) / 3._dp
   g_fac(4) = 8._dp * g2Mi2(1) / 15._dp + 8._dp * g2Mi2(3) / 3._dp
   g_fac(5) = 1._dp * g2Mi2(1) / 30._dp + 1.5_dp * g2Mi2(2) &
          & + 8._dp * g2Mi2(3) / 3._dp
   g_fac = 4._dp * g_fac
   Do i1=1,3
    DME2(i1,i1) = DME2(i1,i1) - g_fac(1)
    DML2(i1,i1) = DML2(i1,i1) - g_fac(2)
    DMD2(i1,i1) = DMD2(i1,i1) - g_fac(3)
    DMU2(i1,i1) = DMU2(i1,i1) - g_fac(4)
    DMQ2(i1,i1) = DMQ2(i1,i1) - g_fac(5)
   End Do
   TraceA(1) = CTrace(AlaAl)
   TraceA(2) = CTrace(AdaAd)
   TraceA(3) = CTrace(AuaAu)
   DMH2(1) = 3._dp * TraceA(2) + TraceA(1) - g_fac(2)   
   DMH2(2) = 3._dp * TraceA(3) - g_fac(2)   

   If (Nu_R) Then
    DMR2 = 2._dp * aAnuAnu
    DML2 = DML2 + AnuaAnu
    TraceA(4) = CTrace(AnuaAnu)
    DMH2(2) = DMH2(2) + TraceA(4)
   End If

   If (Phi) Then
    DMphi2 = 2._dp * Ah02 + 0.5_dp * Alam2
    DMH2 = DMH2 + Ah02
   End If

   If (S_f) Then
    TraceAY(5) = CTrace(A_haA_h)
    DMphi2 = DMphi2 + TraceA(5)
    DMR2 = DMR2 + A_hA_ha
    DMS2 = A_haA_h
   End If
  End If

  !-------------------------------------
  ! two loop order
  !-------------------------------------
  If (loop_ord.Ge.2) Then
   sumE = ZeroC
   sumL = ZeroC
   sumD = ZeroC
   sumU = ZeroC
   sumQ = ZeroC
   sumH = 0._dp

   g4 = g2**2
   
   YlaYlYlaYl = MatMul2(YlaYl,YlaYl,OnlyDiagonal)
   aYlYlaYlYl = MatMul2(aYlYl,aYlYl,OnlyDiagonal)
   YdaYdYdaYd = MatMul2(YdaYd,YdaYd,OnlyDiagonal)
   aYdYdaYdYd = MatMul2(aYdYd,aYdYd,OnlyDiagonal)
   YuaYuYuaYu = MatMul2(YuaYu,YuaYu,OnlyDiagonal)
   aYuYuaYuYu = MatMul2(aYuYu,aYuYu,OnlyDiagonal)

   sumE = - 2._dp * aYlYlaYlYl &
      & - (6._dp*TraceY(2)+2._dp*TraceY(1) - 6._dp*g2(2) + 1.2_dp*g2(1)) * aYlYl
   sumL = - 2._dp * YlaYlYlaYl &
      & - (3._dp*TraceY(2)+TraceY(1) - 1.2_dp*g2(1) ) * YlaYl
   sumD = - 2._dp * aYdYdaYdYd - 2._dp * MatMul3(Y_da,YuaYu,Y_d,OnlyDiagonal) &
      & - (6._dp*TraceY(2)+2._dp*TraceY(1) - 0.4_dp*g2(1) - 6._dp*g2(2)) * aYdYd
   sumU = - 2._dp * aYuYuaYuYu - 2._dp * MatMul3(Y_ua,YdaYd,Y_u,OnlyDiagonal) &
      & - (6._dp*TraceY(3) + 0.4_dp*g2(1) - 6._dp*g2(2)) * aYuYu
   sumQ = - 2._dp * (YuaYuYuaYu+YdaYdYdaYd)                       & 
      & - (3._dp*TraceY(2) + TraceY(1) - 0.4_dp*g2(1) ) * YdaYd   &
      & - (3._dp*TraceY(3) - 0.8_dp*g2(1) ) * YuaYu

   g_fac(1) = (1.2_dp*b_1(1) + 1.44_dp) * g4(1)
   g_fac(2) = (0.3_dp*b_1(1) + 0.09_dp) * g4(1) + 0.9_dp * g2(1) * g2(2)  &
          & + (1.5_dp*b_1(2) + 2.25_dp) * g4(2)
   g_fac(3) = (2._dp * b_1(1) / 15._dp + 4._dp /225._dp) * g4(1) &
          & + 32._dp / 45._dp * g2(1) * g2(3)                    &
          & + (8._dp * b_1(3) / 3._dp + 64._dp/9._dp) * g4(3)
   g_fac(4) = (8._dp * b_1(1) / 15._dp + 64._dp /225._dp) * g4(1) &
          & + 128._dp / 45._dp * g2(1) * g2(3)                    &
          & + (8._dp * b_1(3) / 3._dp + 64._dp/9._dp) * g4(3)
   g_fac(5) = (1._dp * b_1(1) / 30._dp + 1._dp / 900._dp) * g4(1)   &
          & + (1.5_dp*b_1(2) + 2.25_dp) * g4(2)                     &
          & + 0.1_dp * g2(1) * g2(2) + 8._dp * g2(2) * g2(3)        &
          & + 8._dp / 45._dp * g2(1) * g2(3)                        &
          & + (8._dp * b_1(3) / 3._dp + 64._dp/9._dp) * g4(3)
   Do i1=1,3
    sumE(i1,i1) = sumE(i1,i1) + g_fac(1)
    sumL(i1,i1) = sumL(i1,i1) + g_fac(2)
    sumD(i1,i1) = sumD(i1,i1) + g_fac(3)
    sumU(i1,i1) = sumU(i1,i1) + g_fac(4)
    sumQ(i1,i1) = sumQ(i1,i1) + g_fac(5)
   End Do
   
   TraceY2(1) = cTrace(aYlYlaYlYl)
   TraceY2(2) = cTrace(aYdYdaYdYd)
   TraceY2(3) = cTrace(aYuYuaYuYu)
   TraceY2(4) = cTrace(MatMul2(YdaYd,YuaYu,OnlyDiagonal))

   sumH(1) = - 9._dp * TraceY2(2) - 3._dp * TraceY2(4) - 3._dp * TraceY2(1) &
         & + (-0.4_dp*g2(1) + 16._dp * g2(3)) * TraceY(2)                   &
         & + 1.2_dp * g2(1) * TraceY(1)
   sumH(2) = - 9._dp * TraceY2(3) - 3._dp * TraceY2(4)         &
         & + (0.8_dp*g2(1) + 16._dp * g2(3)) * TraceY(3)
   sumH = sumH + g_fac(2)

   If (Nu_R) Then
    sumR = ZeroC
    YnuaYnuYnuaYnu = MatMul2(YnuaYnu,YnuaYnu,OnlyDiagonal)
    aYnuYnuaYnuYnu = MatMul2(aYnuYnu,aYnuYnu,OnlyDiagonal)

    SumE = SumE - 2._dp * MatMul3(Y_la,YnuaYnu,Y_l,OnlyDiagonal)
    SumL = SumL - 2._dp * YnuaYnuYnuaYnu - (3._dp*TraceY(3)+TraceY(4)) * YnuaYnu
    SumU = SumU - 2._dp * TraceY(4) * aYuYu
    SumQ = SumQ - TraceY(4) * YuaYu
    SumR = - 2._dp * aYnuYnuaYnuYnu                           &
         & - 2._dp * MatMul3(Y_nua,YlaYl,Y_nu,OnlyDiagonal)   &
         & - (6._dp*TraceY(3) + 2._dp*TraceY(4) - 1.2_dp*g2(1) - 6._dp*g2(2)) &
         &    * aYnuYnu

    TraceY2(5) = cTrace(aYnuYnuaYnuYnu)
    TraceY2(6) = cTrace(MatMul2(YlaYl,YnuaYnu,OnlyDiagonal))
    SumH(1) = SumH(1) - TraceY2(6)
    SumH(2) = SumH(2) - 3._dp * TraceY2(5) - TraceY2(6) 
   End If

   GammaE = GammaE + oo16pi2 * sumE
   GammaL = GammaL + oo16pi2 * sumL
   GammaD = GammaD + oo16pi2 * sumD
   GammaU = GammaU + oo16pi2 * sumU
   GammaQ = GammaQ + oo16pi2 * sumQ
   GammaH = GammaH + oo16pi2 * sumH

   If (Nu_R) GammaR = GammaR + oo16pi2 * sumR   

  End If
  !-------------------------------------
  ! final normalization
  !-------------------------------------
  beta = oo16pi2 * beta
  GammaE = oo16pi2 * GammaE
  GammaL = oo16pi2 * GammaL
  GammaD = oo16pi2 * GammaD
  GammaU = oo16pi2 * GammaU
  GammaQ = oo16pi2 * GammaQ
  GammaH = oo16pi2 * GammaH

  GammaAE = oo16pi2 * GammaAE
  GammaAL = oo16pi2 * GammaAL
  GammaAD = oo16pi2 * GammaAD
  GammaAU = oo16pi2 * GammaAU
  GammaAQ = oo16pi2 * GammaAQ
  GammaAH = oo16pi2 * GammaAH

  DME2 = oo16pi2 * GammaAE
  DML2 = oo16pi2 * DML2
  DMD2 = oo16pi2 * DMD2
  DMU2 = oo16pi2 * DMU2
  DMQ2 = oo16pi2 * DMQ2
  DMH2 = oo16pi2 * DMH2

  If (Nu_R) Then
   GammaR = oo16pi2 * GammaR
   GammaAR = oo16pi2 * GammaAR
   DMR2 = oo16pi2 * DMR2
  End If
  If (Phi) Then
   GammaPhi = oo16pi2 * GammaPhi
   GammaAPhi = oo16pi2 * GammaAPhi
   DMphi2 = oo16pi2 * DMphi2
  End If
  If (S_f) Then
   GammaS = oo16pi2 * GammaS
   GammaAS = oo16pi2 * GammaAS
   DMS2 = oo16pi2 * DMS2
  End If

  Iname = Iname - 1

 End Subroutine CalculateAnomalousDimensions


 Subroutine RGEs_func(len,T,GY,F)
 Implicit None
  Integer, Intent(in) :: len
  Real(dp), Intent(in) :: T, GY(len)
  Real(dp), Intent(out) :: F(len)

  Integer :: i1, i2, i_zaehl
  Real(dp) :: Q, GammaH(2), GammaAH(2), MH2(2), DMH2(2), DMphi2, Mphi2
  Real(dp), Dimension(3) :: gauge, beta , Dgauge
  Complex(dp) :: h0, lambda, mu, Dmu, Dh0, Dlam, Ah0, Alam, DAh0, DAlam &
      & , B, DB, GammaPhi, GammaAPhi
  Complex(dp), Dimension(3) :: GammaLH, eps, Deps, betaM, DMi, Mi
  Complex(dp), Dimension(3,3) :: Y_l, Y_d, Y_u, Y_nu, h_pnus, GammaE     &
                  & , GammaR , GammaL, GammaD, GammaU, GammaQ, GammaS  &
                  & , A_l, A_nu, A_d, A_u, Ah, GammaAE     &
                  & , GammaAR , GammaAL, GammaAD, GammaAU, GammaAQ, GammaAS
  Complex(dp), Dimension(3,3) :: DY_l, DY_d, DY_u, DY_nu, DA_l, DA_nu, DA_d  &
     & , DA_u, DA_h, DME2, DMR2, DML2, DMD2, DMU2, DMQ2, DMS2, ME2, MR2, ML2 &
     & , MD2, MU2, MQ2, MS2
  Complex(dp), Dimension(3,3,3) :: LambdaE, LambdaP, LambdaPP

  Iname = Iname + 1
  NameOfUnit(Iname) = "RGEs_func"

  q = t

  OnlyDiagonal = .Not.GenerationMixing
  !-----------------------------------------------------
  ! MSSM gauge and Yukawa couplings
  !-----------------------------------------------------
  gauge = gy(1:3)
  i_zaehl = 4
  Do i1=1,3
   Do i2=1,3
    Y_l(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End Do
  End Do
  If (Nu_R) Then
   Do i1=1,3
    Do i2=1,3
     Y_nu(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
  End If
  Do i1=1,3
   Do i2=1,3
    Y_d(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End Do
  End Do
  Do i1=1,3
   Do i2=1,3
    Y_u(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End Do
  End Do
  If (phi) Then
   h0 = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
   i_zaehl = i_zaehl + 2   
   lambda = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
   i_zaehl = i_zaehl + 2   
  End If
  !----------------------------------
  ! bilinear terms in superpotential 
  !----------------------------------
  If (l_mu) Then
   mu = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
   i_zaehl = i_zaehl + 2
  End If
  If (l_mu.And.LH) Then
   Do i1=1,3
    eps(i1) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End Do
  End If
  !----------------------------------
  ! gauginos
  !----------------------------------
  If (l_lin) Then
   Do i1=1,3
    Mi(i1) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End Do
  End If
  !----------------------------------
  ! A-parameters
  !----------------------------------
  If (l_lin) Then
   Do i1=1,3
    Do i2=1,3
     A_l(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   If (Nu_R) Then
    Do i1=1,3
     Do i2=1,3
      A_l(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
   Do i1=1,3
    Do i2=1,3
     A_d(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     A_u(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   If (Phi) Then
    Ah0 = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
    Alam = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
    i_zaehl = i_zaehl + 2
   End If
   If (S_f) Then
    Do i1=1,3
     Do i2=1,3
      Ah(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
  End If
  If (l_qua) Then
  !----------------------------------
  ! B-parameter 
  !----------------------------------
  B = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
  i_zaehl = i_zaehl + 2
  !----------------------------------
  ! soft SUSY mass squared 
  !----------------------------------
   Do i1=1,3
    Do i2=1,3
     ME2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   If (Nu_R) Then
    Do i1=1,3
     Do i2=1,3
      MR2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
   Do i1=1,3
    Do i2=1,3
     ML2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     MD2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     MU2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     MQ2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   MH2 = gy(i_zaehl:i_Zaehl + 1)
   i_zaehl = i_zaehl + 2
   If (Phi) Then
    Mphi2 = gy(i_zaehl)
    i_zaehl = i_zaehl + 1
   End If
   If (S_f) Then
    Do i1=1,3
     Do i2=1,3
      MS2(i1,i2) = Cmplx( gy(i_zaehl), gy(i_Zaehl + 1),dp )
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
  End If

  If ((I_zaehl-1).Ne.len) Then
   Write(ErrCan,*) "Problem in "//NameOfUnit(Iname)
   Write(ErrCan,*) "GY -> Par, i_zaehl, len",i_zaehl,len
   Call TerminateProgram
  End If

  Call CalculateAnomalousDimensions(gauge**2, Y_l, Y_d, Y_u, Y_nu, h0, lambda &
             & , h_pnus, LambdaE, LambdaP, LambdaPP  &
             & , Mi, A_l, A_nu, A_d, A_u, Ah0, Alam, Ah &
             & , ME2, MR2, ML2, MD2, MU2, MQ2, MH2, Mphi2, MS2 &
             & , beta, GammaE, GammaR, GammaL, GammaD, GammaU, GammaQ, GammaH &
             & , GammaPhi, GammaS, GammaLH  &
             & , betaM, GammaAE, GammaAR, GammaAL, GammaAD &
             & , GammaAU, GammaAQ, GammaAH, GammaAPhi, GammaAS &
             & , DME2, DMR2, DML2, DMD2, DMU2, DMQ2, DMH2, DMphi2, DMS2)
  !-----------------------------
  ! gauge couplings
  !-----------------------------
  Dgauge = gauge * beta
  !-----------------------------
  ! Yukawa couplings couplings
  !-----------------------------
  DY_l = GammaH(1) * Y_l + MatMul2(GammaL,Y_l,OnlyDiagonal) &
     & + MatMul2(Y_l, GammaE, OnlyDiagonal )  
  DY_d = GammaH(1) * Y_d + MatMul2(GammaQ,Y_d,OnlyDiagonal) &
     & + MatMul2(Y_d, GammaD, OnlyDiagonal )
  DY_u = GammaH(2) * Y_u + MatMul2(GammaQ,Y_u,OnlyDiagonal) &
     & + MatMul2(Y_u, GammaU, OnlyDiagonal )
  If (Nu_R) DY_nu = GammaH(2) * Y_nu + MatMul2(GammaL,Y_nu,OnlyDiagonal) &
                & + MatMul2(Y_nu, GammaR, OnlyDiagonal )

  If (Phi) Then
   Dh0 = h0 * (Sum(GammaH) + GammaPhi)
   Dlam = 3._dp * lambda * GammaPhi
  End If
  !---------------------------------
  ! bilinear terms in superpotential
  !---------------------------------
  If (l_mu.And.LH) Then
   Dmu = mu * Sum(GammaH)
   Deps = mu * GammaLH + eps * GammaH(2)
   Do i1=1,3
    Dmu = Dmu + eps(i1) * GammaLH(i1)
    Do i2=1,3
     Deps(i1) = Deps(i1) + eps(i2) * GammaL(i2,i1)
    End Do
   End Do
  Else If (l_mu) Then
   Dmu = mu * Sum(GammaH)
  End If

  !-----------------------------
  ! gauginos
  !-----------------------------
  If (l_lin) DMi = betaM
  !-----------------------------
  ! A-parameters
  !-----------------------------
  If (l_lin) Then
   DA_l = GammaH(1) * A_l + MatMul2(GammaL,A_l,OnlyDiagonal) &
      & + MatMul2(A_l, GammaE, OnlyDiagonal )                &
      & - 2._dp * GammaAH(1) * Y_l                           &
      & - 2._dp * MatMul2(GammaAL,Y_l,OnlyDiagonal)          &
      & - 2._dp * MatMul2(Y_l, GammaAE, OnlyDiagonal )
   DA_d = GammaH(1) * A_d + MatMul2(GammaQ,A_d,OnlyDiagonal) &
      & + MatMul2(A_d, GammaD, OnlyDiagonal )                &
      & - 2._dp * GammaAH(1) * Y_d                           &
      & - 2._dp * MatMul2(GammaAQ,Y_d,OnlyDiagonal)          &
      & - 2._dp * MatMul2(Y_d, GammaAD, OnlyDiagonal )
   DA_u = GammaH(2) * A_u + MatMul2(GammaQ,A_u,OnlyDiagonal) &
      & + MatMul2(A_u, GammaU, OnlyDiagonal )                &
      & - 2._dp * GammaAH(2) * Y_u                           &
      & - 2._dp * MatMul2(GammaAQ,Y_u,OnlyDiagonal)          &
      & - 2._dp * MatMul2(Y_u, GammaAU, OnlyDiagonal )
   If (Nu_R) DA_nu = GammaH(2) * A_nu + MatMul2(GammaL,A_nu,OnlyDiagonal) &
                 & + MatMul2(A_nu, GammaR, OnlyDiagonal )                 &
                 & - 2._dp * GammaAH(2) * Y_nu                            &
                 & - 2._dp * MatMul2(GammaAL,Y_nu,OnlyDiagonal)           &
                 & - 2._dp * MatMul2(Y_nu, GammaAR, OnlyDiagonal )

   If (Phi) Then
    DAh0 = Ah0 * (Sum(GammaH) + GammaPhi)  &
       & - 2._dp * h0 * (Sum(GammaAH) + GammaAPhi)
    DAlam = 3._dp * Alam * GammaPhi - 6._dp * lambda * GammaAPhi
   End If
  End If

  If (l_qua) Then
  !----------------------------------
  ! B-parameter 
  !----------------------------------
   DB = B * mu * Sum(GammaH) - 2._dp * mu * Sum(GammaAH) 
  End If
  !---------------------------------
  ! transfer information to f vector
  !---------------------------------
  f(1:3) = Dgauge
  I_zaehl = 4
  Do i1=1,3
   Do i2=1,3
    f(i_zaehl) = Real(DY_l(i1,i2),dp)
    f(i_zaehl + 1) = Aimag(DY_l(i1,i2))
    i_zaehl = i_zaehl + 2
   End Do
  End Do
  If (Nu_R) Then
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DY_nu(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DY_nu(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
  End Do
  End If
  Do i1=1,3
   Do i2=1,3
    f(i_zaehl) = Real(DY_d(i1,i2),dp)
    f(i_zaehl + 1) = Aimag(DY_d(i1,i2))
    i_zaehl = i_zaehl + 2
   End Do
  End Do
  Do i1=1,3
   Do i2=1,3
    f(i_zaehl) = Real(DY_u(i1,i2),dp)
    f(i_zaehl + 1) = Aimag(DY_u(i1,i2))
    i_zaehl = i_zaehl + 2
   End Do
  End Do

  If (phi) Then
   f(i_zaehl) = Real(Dh0,dp)
   f(i_zaehl + 1) = Aimag(Dh0)
   i_zaehl = i_zaehl + 2
   f(i_zaehl) = Real(Dlam,dp)
   f(i_zaehl + 1) = Aimag(Dlam)
   i_zaehl = i_zaehl + 2
  End If

  !--------------------------------
  ! superpotentail mass parameters
  !--------------------------------
  If (l_mu) Then
   f(i_zaehl) = Real(Dmu,dp)
   f(i_zaehl + 1) = Aimag(Dmu)
   i_zaehl = i_zaehl + 2
  End If
  If (l_mu.And.LH) Then
   Do i1=1,3
    f(i_zaehl) = Real(Deps(i1),dp)
    f(i_zaehl + 1) = Aimag(Deps(i1))
    i_zaehl = i_zaehl + 2
   End Do
  End If
  !----------------------------------
  ! gauginos
  !----------------------------------
  If (l_lin) Then
   Do i1=1,3
    f(i_zaehl) = Real(DMi(i1),dp)
    f(i_zaehl + 1) = Aimag(DMi(i1))
    i_zaehl = i_zaehl + 2
   End Do
  End If
  !----------------------------------
  ! A-parameters
  !----------------------------------
  If (l_lin) Then
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DA_l(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DA_l(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   If (Nu_R) Then
    Do i1=1,3
     Do i2=1,3
      f(i_zaehl) = Real(DA_nu(i1,i2),dp)
      f(i_zaehl + 1) = Aimag(DA_nu(i1,i2))
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DA_d(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DA_d(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DA_u(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DA_u(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do

   If (phi) Then
    f(i_zaehl) = Real(DAh0,dp)
    f(i_zaehl + 1) = Aimag(DAh0)
    i_zaehl = i_zaehl + 2
    f(i_zaehl) = Real(DAlam,dp)
    f(i_zaehl + 1) = Aimag(DAlam)
    i_zaehl = i_zaehl + 2
   End If
   If (S_f) Then
    Do i1=1,3
     Do i2=1,3
      f(i_zaehl) = Real(DA_h(i1,i2),dp)
      f(i_zaehl + 1) = Aimag(DA_h(i1,i2))
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
  End If

  If (l_qua) Then
  !----------------------------------
  ! B-parameter 
  !----------------------------------
   f(i_zaehl) = Real(DB,dp)
   f(i_zaehl + 1) = Aimag(DB)
   i_zaehl = i_zaehl + 2
  !----------------------------------
  ! soft SUSY mass squared 
  !----------------------------------
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DME2(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DME2(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   If (Nu_R) Then
    Do i1=1,3
     Do i2=1,3
      f(i_zaehl) = Real(DMR2(i1,i2),dp)
      f(i_zaehl + 1) = Aimag(DMR2(i1,i2))
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DML2(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DML2(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DMD2(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DMD2(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DMU2(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DMU2(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   Do i1=1,3
    Do i2=1,3
     f(i_zaehl) = Real(DMQ2(i1,i2),dp)
     f(i_zaehl + 1) = Aimag(DMQ2(i1,i2))
     i_zaehl = i_zaehl + 2
    End Do
   End Do
   f(i_zaehl:i_Zaehl + 1) = DMH2
   i_zaehl = i_zaehl + 2

   If (phi) Then
    f(i_zaehl) = DMphi2
    i_zaehl = i_zaehl + 1
   End If
   If (S_f) Then
    Do i1=1,3
     Do i2=1,3
      f(i_zaehl) = Real(DMS2(i1,i2),dp)
      f(i_zaehl + 1) = Aimag(DMS2(i1,i2))
      i_zaehl = i_zaehl + 2
     End Do
    End Do
   End If
  End If

  If ((I_zaehl-1).Ne.len) Then
   Write(ErrCan,*) "Problem in "//NameOfUnit(Iname)
   Write(ErrCan,*) "Par -> F, i_zaehl, len",i_zaehl,len
   Call TerminateProgram
  End If

  Iname = Iname - 1

 End Subroutine RGEs_func


 Subroutine SetModel(i_loop, i_model, Wmass, Vmass, Vmass2, rp_bar)
 Implicit None
  Integer, Intent(in) :: i_model     ! model
  Integer, Intent(in) :: i_loop      ! loop order
  Logical, Intent(in) :: Wmass, Vmass, Vmass2
  Logical, Intent(in), Optional :: rp_bar
 
  Iname = Iname + 1
  NameOfUnit(Iname) = "SetModel"
  !------------------------------
  ! initialisation
  !------------------------------
  Nu_R = .False.  ! include contributions due to nu_R
  Phi = .False.   ! include contributions due to phi field
  S_f = .False.   ! include contributions due to S field
      ! include contributions due to explicit lepton number violationg terms
  LH = .False.
   ! include contributions due to explicit baryon number violationg terms
  UDD = .False.
  l_mu = Wmass    ! calculate bilinear terms of the superpotential
  l_lin = Vmass   ! calculate dim [m] terms of the SOFT SUSY potential
  l_qua = Vmass2   ! calculate dim [m^2] terms of the SOFT SUSY potential
  !----------------------------------
  ! loop order
  !----------------------------------
  loop_ord = i_loop  
  If (loop_ord.Gt.1) Then
    Write(*,*)  "only 1-loop RGEs are included up to now"
    Write(ErrCan,*) "only 1-loop RGEs are included up to now"
  End If

  loop_ord = 1

  Select Case(i_model)
  Case(1) ! MSSM
  Case(2) ! MSSM + nu_R
   Nu_R = .True.
  Case(3) ! bilinear RP model
    LH = .True.
  Case(4) ! bilinear + trilinear RP model
   LH = .True.
   If (Present(rp_bar)) UDD = rp_bar 
  Case(5) ! NMSSM
   Phi = .True. 
  Case(6) ! sponaneous RP model
   Phi = .True. 
   Nu_R = .True.
   S_f = .True.
  Case default
   Write(*,*) "model",i_model,"is not defined"
   Write(ErrCan,*) "model",i_model,"is not defined"
   Call TerminateProgram()
  End Select

  Iname = Iname - 1

 End Subroutine SetModel


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
  sumT1 = aYeYe + 6._dp * aYTYT + aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)                 &
        & + Matmul2(Transpose(aYeYe+aYZYZ),YT,OnlyDiagonal)

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
            &        - 0.6_dp * gauge2(1) - 8._dp * gauge2(2)  )
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
     &    * Real(cTrace(Mu)) / 3._dp                           &
     & + (8._dp * gauge2(3) + 0.4_dp*gauge2(1) )               &
     &    * Real(cTrace(Md)) / 3._dp                           &
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
   S1 = mH(2) - mH(1)
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

   betaMT2(1) = MT(1) * (lam12 + TraceY(2)) + 2._dp * mH(1) * lam12   &
            & + 2._dp * Real( cTrace(aYTMlYT),dp ) + TraceA(2) + Alam12 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2(2) = MT(2) * lam22 + 2._dp * mH(2) * lam22 + Alam22 &
            & - 2.4_dp * AbsGM2(1) - 8._dp * AbsGM2(2)

   betaMT2 = 2._dp * betaMT2
 

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
  sumT1 = aYeYe + 6._dp * aYTYT + aYZYZ
  Do i1=1,3
   sumT1(i1,i1) = sumT1(i1,i1) + diagonal(2,1)
  End Do

  betaYT1 = Matmul2(YT,sumT1,OnlyDiagonal)                 &
        & + Matmul2(Transpose(aYeYe+aYZYZ),YT,OnlyDiagonal)

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
  sumT1 = aYeYe + 9._dp * aYTYT + 3._dp * aYZYZ
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
  sumZ1 = sumZ1 + 2._dp * aYZYZ
  betaAZ1 = MatMul2(AZ,sumZ1,OnlyDiagonal)                                   &
    & + 2._dp * MatMul2(2._dp*YSaYs + YdaYd + 2._dp*YZaYZ, AZ,OnlyDiagonal)  &
    & + 6._dp * MatMul2(YZ, aYTAT, OnlyDiagonal)

  diagonal(5,1) = - 2._dp * ( c1_1(2,1) * g2Mi(1) + c1_1(2,2) * g2Mi(2) &
                &           + c1_1(2,3) * g2Mi(3) )

  betaAZ1 = betaAZ1 + diagonal(5,1) * YZ                         &
        & + 2._dp * MatMul2(AdaYd + 4._dp * ASaYS, YZ,OnlyDiagonal)

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
   S1 = mH(2) - mH(1)
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

   betaMZ2(1) = 2._dp * MZ(1) * TraceY(5) + TraceA(5) + betaMZ2(2) &
            & + Real( cTrace(aYZMdYZ) + cTrace(YZMlaYZ) ,dp ) 

   betaMS2(2) = - (3.2_dp * AbsGM2(1) + 40._dp * AbsGM2(3) ) / 3._dp

   betaMS2(1) = MS(1) * TraceY(6) + 2._dp * Real( cTrace(aYSMdYS),dp ) &
            & + 2._dp * TraceA(6) + betaMS2(2)

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
            &        - 0.6_dp * gauge2(1) - 8._dp * gauge2(2)  )
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
   Dgauge = oo16pi2 * gauge * gauge2  &
        & * ( b_1a + oo16pi2 * (Matmul(b_2a,gauge2) - Matmul(a_2a,TraceY(1:4)) ) )
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


End Module RGEs

