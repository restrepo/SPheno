Module DecayFunctions
! comments
!This module is a collection of subroutines with general formulas for the
!calculation of 2-body decays.

! load modules
Use Control
! load modules

Contains


  Subroutine FermionToFermionScalar(mF,mF1,mS,kL,kR,width)
  !-----------------------------------------------------------------------
  ! FermionToFermionScalar(mF,mF1,mS,kL,kR,width) calculates the two body decay
  ! width of fermion decaying to another fermion and a scalar. mF is the mass
  ! of the decaying fermion, mF1 (mS) the mass of the fermion (scalar) in the
  ! final state. kL and kR are the left and right couplings respectively. All
  ! couplings are complex.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
  real(dp), intent(in) :: mF,mF1,mS
  real(dp), intent(out) :: width
  complex(dp), intent(in) :: kL,kR

  real(dp) :: mFsq,mF1sq,mSsq,kappa

  if ( abs(mF).le.( abs(mF1) + mS ) ) then
   width = 0._dp

  elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
   width = 0._dp

  elseif ((mF1.eq.0._dp).and.(mS.eq.0._dp)) then
   If (kL.eq.0._dp) Then
    width = Abs(kR)**2 * mF * oo32Pi
   Else If (kR.eq.0._dp) Then
    width = Abs(kL)**2 * mF * oo32Pi
   Else
    width = (Abs(kL)**2 + Abs(kR)**2) * mF * oo32Pi
   End If

  elseif (mF1.eq.0._dp) then
   mFsq = mF * mF
   mSsq = mS * mS
   If (kL.eq.0._dp) Then
    width =  oo32pi * Abs(kR)**2 * (mFsq - mSsq)**2 /Abs(mF)**3 
   Else If (kR.eq.0._dp) Then
    width =  oo32pi * Abs(kL)**2 * (mFsq - mSsq)**2 /Abs(mF)**3 
   Else
    width =  oo32pi * (Abs(kL)**2 + Abs(kR)**2) * (mFsq - mSsq)**2 /Abs(mF)**3 
   End If

  elseif (mS.eq.0._dp) then
   mFsq = mF * mF
   mF1sq = mF1 * mF1
   If (kL.eq.0._dp) Then
    width = oo32Pi * (mFsq - mF1sq) * Abs(kR)**2 * (mFsq + mF1sq) / Abs(mF)**3
   Else If (kR.eq.0._dp) Then
    width = oo32Pi * (mFsq - mF1sq) * Abs(kL)**2 * (mFsq + mF1sq) / Abs(mF)**3
   Else
    width = oo32Pi * (mFsq - mF1sq)                         &
         & * ( (Abs(kL)**2 + Abs(kR)**2) * (mFsq + mF1sq)  &
         &   + 4._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 ) / Abs(mF)**3  
   End If

  else
   mFsq = mF * mF
   mF1sq = mF1 * mF1
   mSsq = mS * mS
   kappa = Sqrt( (mFsq-mF1sq-mSsq)**2 - 4._dp * mF1sq*mSsq )
   If (kL.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kR)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else If (kR.eq.0._dp) Then
    width = oo32Pi * kappa * Abs(kL)**2 * (mFsq + mF1sq - mSsq) / Abs(mF)**3  
   Else
    width = oo32Pi * kappa                                       &
         & * ( (Abs(kL)**2 + Abs(kR)**2) * (mFsq + mF1sq - mSsq) &
         &   + 4._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 )    / Abs(mF)**3  
   End If

  endif

  End Subroutine FermionToFermionScalar


  Subroutine FermionToFermionVectorBoson(mF,mF1,mV,kL,kR,width)
  !-----------------------------------------------------------------------
  ! FermionToFermionVectorBoson calculates the two body decay 
  ! width of fermion decaying to another fermion and a vectorboson. mF is the
  ! mass of the decaying fermion, mF1 (mV) the mass of the fermion 
  ! (vector boson) in the final state. kL and kR are the left and right 
  ! couplings respectively. All couplings are complex.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
   real(dp), intent(in) :: mF,mF1,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: kL,kR

   real(dp) :: mFsq,mF1sq,mVsq,kappa

   if ( abs(mF).le.( abs(mF1) + mV ) ) then
    width = 0._dp

   elseif (mV.eq.0._dp) then
    write(ErrCan,*) 'Warning from subroutine FermionToFermionVectorBoson.'
    write(ErrCan,*) ' Vectorboson mass = 0 has occured, setting width to 0!!!' 
    width = 0._dp

   elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
    width = 0._dp

   elseif (mF1.eq.0._dp) then
    mFsq = mF * mF
    mVsq = mV * mV
    width =  oo16pi * (Abs(kL)**2 + Abs(kR)**2) * (mFsq - mVsq)         &
          & * ( 0.5_dp * ( mFsq**2 /  mVsq + mFsq ) - mVsq ) / Abs(mF)**3 

   else
    mFsq = mF * mF
    mF1sq = mF1 * mF1
    mVsq = mV * mV
    kappa = Sqrt( (mFsq-mF1sq-mVsq)**2 - 4._dp * mF1sq*mVsq )
    width = oo16Pi * kappa                                                   &
        & * ( (Abs(kL)**2 + Abs(kR)**2)                                      &
        &     * ( 0.5_dp * ( (mFsq-mF1sq)**2 / mVsq + mFsq + mF1sq) - mVsq )  &
        &   - 6._dp * Real( Conjg(kL) * kR,dp ) * mF * mF1 )  / Abs(mF)**3  

   endif

  End Subroutine FermionToFermionVectorBoson


  Subroutine ScalarToScalarVectorBoson(mS,mS1,mV,coup,width)
  !-----------------------------------------------------------------------
  ! ScalarToScalarVectorBoson calculates the two body decay 
  ! width of a scalar decaying to another scalar and a vectorboson. mS is the
  ! mass of the decaying scalar, mS1 (mV) the mass of the scalar (vectorboson)
  ! in the final state. coup is a complex coupling.
  ! written by Werner Porod, 3.11.1999
  ! 23.10.2000: porting to f90
  !-----------------------------------------------------------------------
  implicit none
   real(dp), intent(in) :: mS,mS1,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: coup

   real(dp) :: mSsq,mS1sq,mVsq,kappa

   if ( abs(mS).le.( abs(mS1) + mV ) ) then
    width = 0._dp

   elseif (mV.eq.0._dp) then
    write(ErrCan,*) 'Warning from subroutine ScalarToScalarVectorBoson.'
    write(ErrCan,*) 'Vectorboson mass = 0 has occured, setting width to 0!!!' 
    width = 0._dp

  elseif (Abs(coup).eq.0._dp) then
   width = 0._dp

   elseif (mS1.eq.0._dp) then
    mSsq = mS * mS
    mVsq = mV * mV
    width = oo16pi * Abs(coup)**2  * (mSsq - mVsq)**3 / ( mS**3 * mVsq ) 

   else
    mSsq = mS * mS
    mS1sq = mS1 * mS1
    mVsq = mV * mV
    kappa = Sqrt( (mSsq-mS1sq-mVsq)**2 - 4._dp * mS1sq*mVsq )
    width = oo16Pi * Abs(coup)**2 * kappa**3 / ( mS**3 * mVsq ) 

   endif

  End Subroutine ScalarToScalarVectorBoson


 Subroutine ScalarToTwoFermions(mS,mF1,mF2,kL,kR,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 fermions. mS is the mass of the decaying scalar and
 ! mF1 (mF2) are the masses of the fermions.
 ! kL and kR are the left and right couplings respectively. All
 ! couplings can be complex.
 ! written by Werner Porod, 3.11.1999
 ! 9.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mF1,mF2
  real(dp), intent(out) ::width
  complex(dp), intent(in) :: kL,kR

  real(dp) :: mSsq,mF1sq,mF2sq,kappa

  if ( abs(mS).le.( abs(mF1) + abs(mF2) ) ) then
   width = 0._dp

  elseif ((Abs(kL).eq.0._dp).and.(Abs(kR).eq.0._dp)) then
   width = 0._dp

  elseif ((mF1.eq.0._dp).and.(mF2.eq.0._dp)) then
   width = (Abs(kL)**2 + Abs(kR)**2) * mS * oo16Pi

  elseif ((mF1.eq.0._dp).or.(mF2.eq.0._dp)) then
   mSsq = mS * mS
   if (mF1.eq.0._dp) mF1sq = mF2 * mF2
   if (mF2.eq.0._dp) mF1sq = mF1 * mF1
   width = (Abs(kL)**2 + Abs(kR)**2) * (mSsq - mF1sq)**2 * oo16pi / mS**3 

  elseif (mF1.eq.mF2) then
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   kappa = Sqrt( mSsq**2 - 4._dp * mSsq * mF1sq )
   width = oo16Pi * kappa                                        &
     &   *  ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - 2._dp * mF1sq)  &
     &      - 4._dp * Real( Conjg(kL) * kR,dp ) * mF1sq )  / mS**3  

  else
   mSsq = mS * mS
   mF1sq = mF1 * mF1
   mF2sq = mF2 * mF2
   kappa = Sqrt( (mSsq-mF1sq-mF2sq)**2 - 4._dp * mF1sq*mF2sq )
   width = oo16Pi * kappa                                            &
     &   * ( (Abs(kL)**2 + Abs(kR)**2) * (mSsq - mF1sq - mF2sq)      &
     &     - 4._dp * Real( Conjg(kL) * kR,dp ) * mF2 * mF1 ) / mS**3  

  endif

 End Subroutine ScalarToTwoFermions


 Subroutine ScalarToTwoScalars(mS,mS1,mS2,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 scalars. mS is the mass of the decaying scalar and
 ! mS1 (mS2) are the masses of the scalars in the final state.
 ! coup is a complex coupling.
 ! written by Werner Porod, 3.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mS1,mS2
  real(dp), intent(out) ::width
  complex(dp), intent(in) :: coup

  real(dp) :: mSsq,mS1sq,mS2sq,kappa

  if ( abs(mS).le.( abs(mS1) + abs(mS2) ) ) then
   width = 0._dp

  elseif (Abs(coup).eq.0._dp) then
   width = 0._dp

  elseif ((mS1.eq.0._dp).and.(mS2.eq.0._dp)) then
   width = Abs(coup)**2 * oo16Pi / mS 

  elseif ((mS1.eq.0._dp).or.(mS2.eq.0._dp)) then
   mSsq = mS * mS
   if (mS1.eq.0._dp) mS1sq = mS2 * mS2
   if (mS2.eq.0._dp) mS1sq = mS1 * mS1
   width = Abs(coup)**2 * (mSsq - mS1sq) * oo16pi / mS**3 

  else
   mSsq = mS * mS
   mS1sq = mS1 * mS1
   mS2sq = mS2 * mS2
   kappa = Sqrt( (mSsq-mS1sq-mS2sq)**2 - 4._dp * mS1sq*mS2sq )
   width = oo16Pi * Abs(coup)**2 * kappa / mS**3 

  endif

  End Subroutine ScalarToTwoScalars


 Subroutine ScalarToTwoVectorbosons(mS,mV,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 Vectorbosons. mS is the mass of the decaying scalar and
 ! mV are the mass of the vectorboson, coup is a real coupling.
 ! written by Werner Porod, 5.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
  implicit none
   real(dp), intent(in) :: mS,mV
   real(dp), intent(out) :: width
   complex(dp), intent(in) :: coup

   real(dp) :: mSsq,mVsq,x,kappa

   if ( abs(mS).le.( 2._dp * mV ) ) then
    width = 0._dp

   elseif (mV.eq.0._dp) then
    write(ErrCan,*) 'Server warning, in subroutine ScalarToTwoVectorbosons'
    write(ErrCan,*) 'm_V = 0, setting width to 0'
    width = 0._dp

   elseif (Abs(coup).eq.0._dp) then
    width = 0._dp

   else
    mSsq = mS * mS
    mVsq = mV**2
    x = mSsq / mVsq
    kappa = Sqrt( 1._dp - 4._dp / x )
    width = oo64Pi * coup**2 * kappa * (x*x - 4._dp*x + 1.2d1) / mS 

   endif

  End Subroutine ScalarToTwoVectorbosons


 Subroutine ScalarToVectorbosonsVR(mS,mV,coup,width)
 !-----------------------------------------------------------------------
 ! ScalarToTwoFermions calculates the two body decay width of a scalar
 ! decaying to 2 Vectorbosons, where on of them is real and the other is
 ! virtuell. mS is the mass of the decaying scalar and
 ! mV are the mass of the vectorboson. coup is an effective real coupling,
 ! which is the scalar-vectorboson coupling times g times a correction
 ! factor for the Z-boson. This correction factor squared reads:
 ! (m_Z/m_W)**2 * (3.5-20.*sinw2/3.+80.*sinw4/9.) * cosW2
 ! written by Werner Porod, 5.11.1999
 ! 23.10.2000: porting to f90
 !-----------------------------------------------------------------------
 implicit none
  real(dp), intent(in) :: mS,mV,coup
  real(dp), intent(out) :: width

  real(dp) :: RHiggs,x,x2,sum1,sum2,sum3

  width = 0._dp

  if ((mS.lt.mV).or.(Abs(coup).eq.0._dp) ) return
  If (mS.Gt.(2._dp*mV-0.1_dp)) return ! in this case a refined calculation
                                      ! is necessary

  x = (mV / mS )**2
  x2 = x*x
  sum1 = 3._dp * (1._dp-8._dp*x+20._dp*x2) / sqrt(4._dp*x-1._dp)   &
     & * acos( (3._dp*x-1._dp) / (2._dp*sqrt(x**3)) )
  sum2 = -0.5_dp * (1._dp-x) * (2._dp-13._dp*x+47._dp*x2) / x
  sum3 = - 1.5_dp * (1._dp-6._dp*x+4._dp*x2) * log(x)

  Rhiggs = sum1 + sum2 + sum3

  width = 3._dp * oo128pi3 * coup**2 * mS * Rhiggs 

 End Subroutine ScalarToVectorbosonsVR


End Module DecayFunctions

