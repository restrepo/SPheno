Module Slepton3BodyDecays

Use Control
Use Mathematics
Use StandardModel, Only: mf_l

Real(dp), Private :: rl2, rt2, rst2, ri2, rj2

Contains


 Subroutine GMSBintegrands(s, erg)
 Implicit None

  Real(dp), Intent(in) :: s(2)
  Real(dp), Intent(out) :: erg(6)

  Integer :: i1
  Real(dp) :: x, x2, fij

  erg = 0._dp

  Do i1=1,2
   x = s(i1)
   x2 = x**2

   fij = Sqrt(x2-4._dp*rl2) * kappa(1-x+rl2,rt2,rst2)                    &
     & / ( (1-x+rl2)**2 * (ri2-1+x-rl2) * (rj2-1+x-rl2) )

   erg(2) = erg(2) + (x - 2 * rl2) * (1 - x + rl2 + rt2 - rst2) * fij
   fij = (1-x+rl2) * fij 
   erg(1) = erg(1) + (x - 2 * rl2) * (1 - x + rl2 + rt2 - rst2) * fij
   erg(3) = erg(3) + (x - 2 * rl2) * fij
   erg(4) = erg(4) + (1 - x + rl2 + rt2 - rst2) * fij
   erg(5) = erg(5) + fij
   erg(6) = erg(6) + (1-x+rl2) * fij
  End Do

 End Subroutine GMSBintegrands


 Subroutine Slepton_to_Slepton_ll(i, j, k, l, mSl, mN, C_L, C_R, eps, gam)
 Implicit None
  Integer, Intent(in) :: i         ! index of decaying slepton
  Integer, Intent(in) :: j         ! final state slepton
  Integer, Intent(in) :: k         ! index of lepton k
  Integer, Intent(in) :: l         ! index of lepton l
  Real(dp), Intent(in) :: mSl(:)   ! slepton masses
  Real(dp), Intent(in) :: mN(4)    ! neutralino masses
  Complex(dp), Intent(in) :: C_L(:,:,:), C_R(:,:,:) ! LR couplings of sleptons
  Real(dp), Intent(in) :: eps      ! required relative precision
  Real(dp), Intent(out) :: gam(2)  ! partial widths: 1 .. same sign sleptons
                                   !                 2 .. opposite sign sleptons


  Integer :: i1, i2
  Real(dp) :: I_aij(6,4,4), smin, smax

  Iname = Iname + 1
  NameOfUnit(Iname) = "Slepton_to_Slepton_ll"

  gam = 0._dp
  !----------------------------------------------
  ! kinematical bound
  !----------------------------------------------
  If (mSl(i).Lt.(mSl(j)+mf_l(k)+mf_l(l)) ) Then
   Iname = Iname - 1
   Return
  Else If (mSl(i).Gt.Abs(mN(1)) ) Then
   Iname = Iname - 1
   Return
  End If

  !--------------------------------------
  ! kinematical functions
  !--------------------------------------
  rl2 = (mf_l(k) / mSl(i))**2
  rt2 = (mf_l(l) / mSl(i))**2
  rst2 = (mSl(j) / mSl(i) )**2
  smin = 2._dp * mf_l(k) / mSl(j)
  smax = 1._dp + rl2 -  rt2 - rst2 - 2._dp * Sqrt(rl2*rt2)
  Do i1=1,4
   ri2 = (mN(i1) / mSl(i) )**2
   Do i2 = i1,4
    rj2 = (mN(i2) / mSl(i) )**2
    Call DgaussInt(GMSBintegrands,6,smin,smax,I_aij(:,i1,i2),eps)
    I_aij(2,i1,i2) = ri2 * rj2 * I_aij(2,i1,i2)
    I_aij(3,i1,i2) = 2._dp * rt2 * rj2 * I_aij(3,i1,i2)
    I_aij(4,i1,i2) = 2._dp * rl2 * rj2 * I_aij(4,i1,i2)
    I_aij(5,i1,i2) = 2._dp * rl2 * rt2 * ri2 * rj2 * I_aij(5,i1,i2)
    I_aij(6,i1,i2) = 2._dp * rl2 * rt2 * I_aij(6,i1,i2)
    If (i1.Ne.i2) I_aij(:,i2,i1) = I_aij(:,i1,i2)
   End Do
  End Do

  Do i1=1,4
   Do i2=1,4
    gam(1) = gam(1)                                                            &
     & + ( Conjg( C_R(i,k,i1) * C_R(j,l,i2) ) * C_R(i,k,i2) * C_R(j,l,i1)      &
     &   + Conjg( C_L(i,k,i1) * C_L(j,l,i2) ) * C_L(i,k,i2) * C_L(j,l,i1))     &
     &  * I_aij(1,i1,i2)                                                       &
     & + ( Conjg( C_R(i,k,i1) * C_L(j,l,i2) ) * C_R(i,k,i2) * C_L(j,l,i1)      &
     &   + Conjg( C_L(i,k,i1) * C_R(j,l,i2) ) * C_L(i,k,i2) * C_R(j,l,i1))     &
     &  * I_aij(2,i1,i2)                                                       &
     & + Real( Conjg( C_R(i,k,i1) * C_R(j,l,i2) ) * C_R(i,k,i2) * C_L(j,l,i1)  &
     &       + Conjg( C_L(i,k,i1) * C_L(j,l,i2) ) * C_L(i,k,i2) * C_R(j,l,i1)) &
     &      * 2._dp * I_aij(3,i1,i2)                                           &
     & - Real( Conjg( C_R(i,k,i1) * C_L(j,l,i2) ) * C_L(i,k,i2) * C_L(j,l,i1)  &
     &       + Conjg( C_L(i,k,i1) * C_R(j,l,i2) ) * C_R(i,k,i2) * C_R(j,l,i1)) &
     &      * 2._dp * I_aij(4,i1,i2)                                           &
     & - Real( Conjg( C_L(i,k,i1) * C_L(j,l,i2) ) * C_R(i,k,i2) * C_R(j,l,i1)) &
     &      * 4._dp * I_aij(5,i1,i2)                                           &
     & - Real( Conjg( C_R(i,k,i1) * C_R(j,l,i2) ) * C_L(i,k,i2) * C_L(j,l,i1)) &
     &      * 4._dp * I_aij(6,i1,i2) 

    gam(2) = gam(2)                                                            &
     & + ( Conjg( C_R(i,k,i1) * C_L(j,l,i2) ) * C_R(i,k,i2) * C_L(j,l,i1)      &
     &   + Conjg( C_L(i,k,i1) * C_R(j,l,i2) ) * C_L(i,k,i2) * C_R(j,l,i1))     &
     &  * I_aij(1,i1,i2)                                                       &
     & + ( Conjg( C_R(i,k,i1) * C_R(j,l,i2) ) * C_R(i,k,i2) * C_R(j,l,i1)      &
     &   + Conjg( C_L(i,k,i1) * C_L(j,l,i2) ) * C_L(i,k,i2) * C_L(j,l,i1))     &
     &  * I_aij(2,i1,i2)                                                       &
     & + Real( Conjg( C_R(i,k,i1) * C_L(j,l,i2) ) * C_R(i,k,i2) * C_R(j,l,i1)  &
     &       + Conjg( C_L(i,k,i1) * C_R(j,l,i2) ) * C_L(i,k,i2) * C_L(j,l,i1)) &
     &      * 2._dp * I_aij(3,i1,i2)                                           &
     & - Real( Conjg( C_R(i,k,i1) * C_R(j,l,i2) ) * C_L(i,k,i2) * C_R(j,l,i1)  &
     &       + Conjg( C_L(i,k,i1) * C_L(j,l,i2) ) * C_R(i,k,i2) * C_L(j,l,i1)) &
     &      * 2._dp * I_aij(4,i1,i2)                                           &
     & - Real( Conjg( C_L(i,k,i1) * C_R(j,l,i2) ) * C_R(i,k,i2) * C_L(j,l,i1)) &
     &      * 4._dp * I_aij(5,i1,i2)                                           &
     & - Real( Conjg( C_R(i,k,i1) * C_L(j,l,i2) ) * C_L(i,k,i2) * C_R(j,l,i1)) &
     &      * 4._dp * I_aij(6,i1,i2) 
   End Do
  End Do

  If (j.Eq.k) gam(2) = 0.5_dp * gam(2)

  gam = oo512pi3 * mSl(i) * gam

  Iname = Iname - 1

 End Subroutine Slepton_to_Slepton_ll


End Module Slepton3bodyDecays
