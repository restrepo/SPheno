Module ThreeBodyPhaseSpaceS
! comments
!-----------------------------------------------------------
! contains the phase space function needed for 3-body decays
! involving fermions and bosons in intial and/or final state 
!-----------------------------------------------------------

! load modules
Use Control
! load modules

Contains

 Complex(dp) Function I2dsmg1tmg2(s, tmin, tmax, m12, g1, m22, g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     1                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
  Implicit None

  Real(dp), Intent(in) :: s, tmin, tmax, m12, m22, g1, g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2dsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0) Then
   ReT = Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   ReT = 0.5_dp * Log(( g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    ImT = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    ImT = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    ImT = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
  End If

  I2dsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

  End Function I2dsmg1tmg2


 Complex(dp) Function I2dtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     1                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3),g12,g22, x1, x2, y1, y2

  I2dtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2dtm12g12 = ( Log(Abs( (tmax-m12) / (tmin-m12) ) )          &
              & - Log(Abs( (tmax-m22) / (tmin-m22) ) ) ) / (m12 - m22)
  Else
   g12 = g1 * g1
   g22 = g2 * g2
   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * (m12-m22)
   fa(3) = (g1+g2)
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2dtm12g12 =                                                            &
     &   fa(2) * ( Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )    &
     &          - Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22) )  )  &
     & - fa(3) * ( y1 + y2  )
   I2dtm12g12 = I2dtm12g12 / fa(1)

  End If

 End Function I2dtm12g12


 Real(dp) Function I2dtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (            1
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin, tmax, m12, g12
  Real(dp) :: g1, x1, x2

  If (g12.Le.0._dp) Then
   I2dtm1g1 = 1._dp / (tmin-m12) - 1._dp / (tmax-m12)
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   I2dtm1g1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) I2dtm1g1 = I2dtm1g1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) I2dtm1g1 = I2dtm1g1 - pi
   I2dtm1g1 = I2dtm1g1 / g1
  End If

  End Function I2dtm1g1


 Complex(dp) Function I2t2dsmg1tmg2(s,tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t^2                 )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s,tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2t2dsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0._dp) Then
   ReT = 0.5_dp * ( tmax*(2._dp*m22+tmax) - tmin*(2._dp*m22+tmin)) &
     & + m22**2 * Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   fa(1) = Log((g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   ReT = 0.5_dp * ( tmax*(2._dp*m22+tmax) - tmin*(2._dp*m22+tmin)) &
     & + 0.5_dp * (m22**2 - g2**2) * fa(1)
   ImT = -g2 * (tmax - tmin + m22 * fa(1))
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    fa(2) = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    fa(2) = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    fa(2) = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
   ReT = ReT + 2._dp * g2 * m22 * fa(2)
   ImT = ImT + (m22**2 - g2**2) * fa(2)
  End If

  I2t2dsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

 End Function I2t2dsmg1tmg2


 Complex(dp) Function I2tdtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(4),g12,g22, x1, x2, y1, y2

  I2tdtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2tdtm12g12 = ( m12 * Log(Abs( (tmax-m12) / (tmin-m12) ) )      &
               & - m22 * Log(Abs( (tmax-m22) / (tmin-m22) ) ) ) / (m12 - m22)
  Else
   g12 = g1 * g1
   g22 = g2 * g2
   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * (m12 * (m12-m22) + g1 *(g1+g2) )
   fa(3) = 0.5_dp * (m22 * (m22-m12) + g2 *(g1+g2) )
   fa(4) = (g1*m22+g2*m12)
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2tdtm12g12 = fa(2) * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) ) &
             & + fa(3) * Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22) ) &
             & - fa(4) * ( y1 + y2  )
   I2tdtm12g12 = I2tdtm12g12 / fa(1)
  End If

 End Function I2tdtm12g12


 Complex(dp) Function I2t2dtm12g12(tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                    t^2                  )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (t - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None

  Real(dp), Intent(in) :: tmin,tmax,m12,m22,g1,g2
  Real(dp) :: m14,m24,fa(5),g12,g22, x1, x2, y1, y2

  I2t2dtm12g12 = ZeroC

  If ((g1.Eq.0._dp).Or.(g2.Eq.0._dp)) Then
   I2t2dtm12g12 = tmax - tmin                                           &
              & + ( m12**2 * Log(Abs( (tmax-m12) / (tmin-m12) ) )       &
              &   - m22**2 * Log(Abs( (tmax-m22) / (tmin-m22) ))) / (m12 - m22)
  Else
   m14 = m12 * m12
   m24 = m22 * m22
   g12 = g1 * g1
   g22 = g2 * g2

   fa(1) = (m12-m22)**2 + (g1+g2)**2
   fa(2) = 0.5_dp * ( (m14-g12)*(m12-m22) + 2._dp* m12* g1* (g1+g2) )
   fa(3) = -( (m14-g12)*(g1+g2) - 2._dp* m12* g1* (m12-m22) )
   fa(4) = 0.5_dp * ( (m24-g22)*(m22-m12) + 2._dp* m22* g2* (g1+g2) )
   fa(5) = -( (m24-g22)*(g1+g2) + 2._dp* m22* g2* (m12-m22) )
   x1 = g1 / (tmax-m12)
   x2 = g1 / (tmin-m12)
   y1 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y1 = y1 - pi
   x1 = g2 / (tmax-m22)
   x2 = g2 / (tmin-m22)
   y2 = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y2 = y2 - pi

   I2t2dtm12g12 = tmax - tmin                                                 &
              & + fa(2) * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12)) &
              & + fa(3) * y1                                                  &
              & + fa(4) * Log( ((tmax-m22)**2 + g22) / ((tmin-m22)**2 + g22)) &
              & + fa(5) * y2
   I2t2dtm12g12 = I2t2dtm12g12 / fa(1)
  End If

 End Function I2t2dtm12g12



 Real(dp) Function I2t2dtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (           t^2
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp),Intent(in) :: tmin, tmax, m12, g12
  Real(dp) :: g1, x1, x2, y

  If (g12.Le.0._dp) Then
   I2t2dtm1g1 = m12**2 * ( 1._dp/(tmin-m12) - 1._dp/(tmax-m12) )    &
            & + 2._dp * m12 * Log(Abs( (tmax-m12) / (tmin-m12) ) )  &
            & + tmax - tmin
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   y = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y - pi

   I2t2dtm1g1 = (m12**2-g12) * y / g1 &
            & + m12 * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )   &
            & + tmax - tmin
  End If

 End Function I2t2dtm1g1


 Complex(dp) Function I2tdsmg1tmg2(s,tmin,tmax,m12,g1,m22,g2)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   ( tmax                                        )
 !   (   (                     t                   )
 ! Re(   | dt -----------------------------------  )
 !   (   )     (s - m12 + I g1) (t - m22 - I g2)   )
 !   ( tmin                                        )
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: s,tmin,tmax,m12,m22,g1,g2
  Real(dp) :: fa(3), ReS, ImS, ReT, ImT

  I2tdsmg1tmg2 = ZeroC

  If (g1.Eq.0._dp) Then
   ReS = 1._dp / (s-m12)
   ImS = 0._dp
  Else
   ReS = (s-m12) / ( (s-m12)**2 + g1**2)
   ImS = - g1 / ( (s-m12)**2 + g1**2)
  End If

  If (g2.Eq.0._dp) Then
   ReT = tmax - tmin + m22 * Log(Abs(tmax - m22) / Abs(tmin - m22) )
   ImT = 0._dp
  Else
   fa(1) = 0.5_dp * Log((g2**2 + (tmax - m22)**2) / ( g2**2 + (tmin - m22)**2))
   ReT = tmax - tmin + m22 * fa(1)
   ImT = -g2 * fa(1)
   fa(1) = (tmin - m22) / g2
   fa(2) = (tmax - m22) / g2
   fa(3) = fa(1) * fa(2)
   If (fa(3).Gt.-1._dp) Then
    fa(2) = Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else If (fa(1).Gt.0._dp) Then
    fa(2) = Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   Else
    fa(2) = -Pi + Atan( (fa(1) - fa(2)) / (1+fa(3)) )
   End If 
   ReT = ReT + g2 * fa(2)
   ImT = ImT + m22 * fa(2)
  End If

  I2tdsmg1tmg2 = Cmplx(ReS * ReT - ImS * ImT , ReT * ImS + ReS * ImT,dp)

 End Function I2tdsmg1tmg2


 Real(dp) Function I2tdtm1g1(tmin, tmax, m12, g12)
 !--------------------------------------------------------------------
 ! Computes the following integral:
 !
 !   tmax
 !     (            t
 !     | dt -----------------
 !     )     (t-m12)^2 + g12
 !   tmin
 !
 ! written by Werner Porod, 9.1.1996
 !--------------------------------------------------------------------
 Implicit None
  Real(dp), Intent(in) :: tmin,tmax,m12,g12
  Real(dp) :: g1, x1, x2, y


  If (g12.Le.0._dp) Then
   I2tdtm1g1 = m12 / (tmin-m12) - m12 / (tmax-m12)   &
           & + Log(Abs( (tmax-m12) / (tmin-m12) ) )
  Else
   g1 = Sqrt(g12)
   x1 = (tmax-m12)/g1
   x2 = (tmin-m12)/g1
   y = Atan((x1-x2)/(1._dp + x1 * x2))
   If ((x1.Gt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y + pi
   If ((x1.Lt.0._dp).And.(x1*x2.Lt.-1._dp)) y = y - pi

   I2tdtm1g1 = m12 * y / g1          &
           & + 0.5_dp * Log( ((tmax-m12)**2 + g12) / ((tmin-m12)**2 + g12) )
  End If

 End Function I2tdtm1g1

End Module ThreeBodyPhaseSpaceS

