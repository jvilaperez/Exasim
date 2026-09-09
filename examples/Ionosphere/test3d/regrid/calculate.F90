module calculate_module

  use prec,only:rp

  implicit none

  real(kind=rp),parameter :: rmass_o1 = 16, rmass_o2 = 32, &
    rmass_n1 = 14, rmass_n2 = 28, rmass_no = 30, rmass_he = 4

  contains
!-----------------------------------------------------------------------
  elemental subroutine collision_frequency(tr,o2,o1,he,n2, &
    nu_op,nu_o2p,nu_np,nu_n2p,nu_nop)

    real(kind=rp),intent(in) :: tr,o2,o1,he,n2
    real(kind=rp),intent(out) :: nu_op,nu_o2p,nu_np,nu_n2p,nu_nop

    real(kind=rp) :: sqrttr,lgtr,op_r,o2p_r,n2p_r, &
      op_nr,o2p_nr,np_nr,n2p_nr,nop_nr

    op_nr  = 6.64_rp*o2              + 1.32_rp*he + 6.82_rp*n2
    o2p_nr =              2.31_rp*o1 + 0.70_rp*he + 4.13_rp*n2
    np_nr  = 7.25_rp*o2 + 4.42_rp*o1 + 1.49_rp*he + 7.47_rp*n2
    n2p_nr = 4.49_rp*o2 + 2.58_rp*o1 + 0.79_rp*he
    nop_nr = 4.27_rp*o2 + 2.44_rp*o1 + 0.74_rp*he + 0.69_rp*n2

    sqrttr = sqrt(tr)
    lgtr = log10(tr)
    op_r  = 3.67_rp*sqrttr*(1 - 0.064_rp*lgtr)**2
    o2p_r = 2.59_rp*sqrttr*(1 - 0.073_rp*lgtr)**2
    n2p_r = 5.14_rp*sqrttr*(1 - 0.069_rp*lgtr)**2

    nu_op = 1e-10_rp*op_nr + 1e-11_rp*op_r*o1
    nu_o2p = 1e-10_rp*o2p_nr + 1e-11_rp*o2p_r*o2
    nu_np = 1e-10_rp*np_nr
    nu_n2p = 1e-10_rp*n2p_nr + 1e-11_rp*n2p_r*n2
    nu_nop = 1e-10_rp*nop_nr

  endsubroutine collision_frequency
!-----------------------------------------------------------------------
  elemental subroutine equivalent_efield( &
    E1,E2,E3,u1,u2,u3,B1,B2,B3,Eeq1,Eeq2,Eeq3)

    real(kind=rp),intent(in) :: E1,E2,E3,u1,u2,u3,B1,B2,B3
    real(kind=rp),intent(out) :: Eeq1,Eeq2,Eeq3

    Eeq1 = E1 + u2*B3 - u3*B2
    Eeq2 = E2 + u3*B1 - u1*B3
    Eeq3 = E3 + u1*B2 - u2*B1

  endsubroutine equivalent_efield
!-----------------------------------------------------------------------
  elemental function ion_heat(ped,Eeq1,Eeq2,Eeq3) result(heat_i)

    real(kind=rp),intent(in) :: ped,Eeq1,Eeq2,Eeq3
    real(kind=rp) :: heat_i

    heat_i = ped * (Eeq1**2 + Eeq2**2 + Eeq3**2)

  endfunction ion_heat
!-----------------------------------------------------------------------
  elemental function electron_heat_coeff(r) result(coeff)

    use util_module,only:isclose

    real(kind=rp),intent(in) :: r
    real(kind=rp) :: coeff

    real(kind=rp),dimension(0:6),parameter :: c = &
      (/5.342_rp, 1.056_rp, -4.392e-2_rp, -5.9e-2_rp, &
        -9.346e-3_rp, -5.755e-4_rp, -1.249e-5_rp/)
    integer :: i
    real(kind=rp) :: lnr,s

    s = c(0)

    if (.not. isclose(r,0.0_rp)) then
      lnr = log(r)
      do i = 1,6
        s = s + c(i)*lnr**i
      enddo
    endif

    coeff = exp(s)

  endfunction electron_heat_coeff
!-----------------------------------------------------------------------
  elemental function transfer_electron_ion(te,ne,op,o2p,np,n2p,nop) result(qei)

    real(kind=rp),intent(in) :: te,ne,op,o2p,np,n2p,nop
    real(kind=rp) :: qei

    real(kind=rp),parameter :: lnlambda = 15

    qei = 5.12e-7_rp*ne/te**1.5_rp*lnlambda* &
      (op/rmass_o1+o2p/rmass_o2+np/rmass_n1+n2p/rmass_n2+nop/rmass_no)

  endfunction transfer_electron_ion
!-----------------------------------------------------------------------
  elemental function transfer_ion_neutral( &
    tr,o2,o1,he,n2,op,o2p,np,n2p,nop) result(qin)

    real(kind=rp),intent(in) :: tr,o2,o1,he,n2,op,o2p,np,n2p,nop
    real(kind=rp) :: qin

    real(kind=rp),parameter :: ao2 = 1.6_rp, ao1 = 0.89_rp, ahe = 0.21_rp, an2 = 1.76_rp
    real(kind=rp) :: nonresonant,resonant

    nonresonant = &
      op *o2*sqrt(ao2*rmass_o1*rmass_o2/(rmass_o1+rmass_o2)**3)+ &
      op *he*sqrt(ahe*rmass_o1*rmass_he/(rmass_o1+rmass_he)**3)+ &
      op *n2*sqrt(an2*rmass_o1*rmass_n2/(rmass_o1+rmass_n2)**3)+ &
      o2p*o1*sqrt(ao1*rmass_o2*rmass_o1/(rmass_o2+rmass_o1)**3)+ &
      o2p*he*sqrt(ahe*rmass_o2*rmass_he/(rmass_o2+rmass_he)**3)+ &
      o2p*n2*sqrt(an2*rmass_o2*rmass_n2/(rmass_o2+rmass_n2)**3)+ &
      np *o2*sqrt(ao2*rmass_n1*rmass_o2/(rmass_n1+rmass_o2)**3)+ &
      np *o1*sqrt(ao1*rmass_n1*rmass_o1/(rmass_n1+rmass_o1)**3)+ &
      np *he*sqrt(ahe*rmass_n1*rmass_he/(rmass_n1+rmass_he)**3)+ &
      np *n2*sqrt(an2*rmass_n1*rmass_n2/(rmass_n1+rmass_n2)**3)+ &
      n2p*o2*sqrt(ao2*rmass_n2*rmass_o2/(rmass_n2+rmass_o2)**3)+ &
      n2p*o1*sqrt(ao1*rmass_n2*rmass_o1/(rmass_n2+rmass_o1)**3)+ &
      n2p*he*sqrt(ahe*rmass_n2*rmass_he/(rmass_n2+rmass_he)**3)+ &
      nop*o2*sqrt(ao2*rmass_no*rmass_o2/(rmass_no+rmass_o2)**3)+ &
      nop*o1*sqrt(ao1*rmass_no*rmass_o1/(rmass_no+rmass_o1)**3)+ &
      nop*he*sqrt(ahe*rmass_no*rmass_he/(rmass_no+rmass_he)**3)+ &
      nop*n2*sqrt(an2*rmass_no*rmass_n2/(rmass_no+rmass_n2)**3)

    resonant = 1.4_rp*o2p*o2 + 2.1_rp*op*o1 + 2.7_rp*n2p*n2

    qin = 6.8e-13_rp*nonresonant + 1e-15_rp*sqrt(2*tr)*resonant

  endfunction transfer_ion_neutral
!-----------------------------------------------------------------------
  elemental function transfer_electron_neutral(tn,o2,o1,he,n2,ne,te) result(qen)

    use util_module,only:isclose

    real(kind=rp),intent(in) :: tn,o2,o1,he,n2,ne,te
    real(kind=rp) :: qen

    real(kind=rp) :: dt,sqrtte,a,r,cool_o2,cool_o1,cool_he,cool_n2

    dt = te-tn
    sqrtte = sqrt(te)

    if (te <= 1000) then
      a = 5.71e-8_rp*exp(-3352.6_rp/te)
    elseif (te < 2000) then
      a = 2e-7_rp*exp(-4605.2_rp/te)
    else
      a = 2.53e-6_rp*sqrtte*exp(-17620/te)
    endif

    if (isclose(dt,0.0_rp)) then
      r = 3200/tn**2
    else
      r = (1-exp(-3200*dt/(te*tn)))/dt
    endif

    cool_o2 = 1.21e-18_rp*(sqrtte + 3.6e-2_rp*te) + &
      6.9e-14_rp/sqrtte + 3.125e-21_rp*te**2
    cool_o1 = 7.9e-19_rp*sqrtte*(1 + 5.7e-4_rp*te) + &
      3.4e-12_rp*max(1 - 7e-5_rp*te,0.0_rp)*(150/te + 0.4_rp)/tn
    cool_he = 2.46e-17_rp*sqrtte
    cool_n2 = 1.77e-19_rp*te*max(1 - 1.21e-4_rp*te,0.0_rp) + &
      2.9e-14_rp/sqrtte + 1.3e-4_rp*a*r
    qen = ne*(cool_o2*o2 + cool_o1*o1 + cool_he*he + cool_n2*n2)

  endfunction transfer_electron_neutral
!-----------------------------------------------------------------------
  elemental function ion_conductivity(op,ne) result(kappa_i)

    use util_module,only:isclose

    real(kind=rp),intent(in) :: op,ne
    real(kind=rp) :: kappa_i

    if (isclose(ne,0.0_rp)) then
      kappa_i = 0
    else
      kappa_i = 1.15e4_rp*op/ne
    endif

  endfunction ion_conductivity
!-----------------------------------------------------------------------
  elemental function electron_conductivity(te,ne,o2,o1,he,n2) result(kappa_e)

    use util_module,only:isclose

    real(kind=rp),intent(in) :: te,ne,o2,o1,he,n2
    real(kind=rp) :: kappa_e

    real(kind=rp) :: sqrtte,q

    sqrtte = sqrt(te)

    q = (2.2_rp + 7.92e-2_rp*sqrtte)*o2 + &
        1.1_rp*(1 + 5.7e-4_rp*te)*o1 + &
        5.6_rp*he + &
        max(0.282_rp - 3.41e-5_rp*te,0.0_rp)*sqrtte*n2

    if (isclose(ne,0.0_rp)) then
      kappa_e = 0
    else
      kappa_e = 7.5e21_rp/(1e16_rp + 3.22e4_rp*te**2*q/ne)
    endif

  endfunction electron_conductivity
!-----------------------------------------------------------------------
  elemental subroutine calculate(tn,ti,te,o2,o1,he,n2, &
    op,o2p,np,n2p,nop,ped,qsum,E1,E2,E3,u1,u2,u3,B1,B2,B3, &
    nu_op,nu_o2p,nu_np,nu_n2p,nu_nop, &
    heat_i,heat_e,qei,qin,qen,kappa_i,kappa_e)

    real(kind=rp),intent(in) :: tn,ti,te,o2,o1,he,n2, &
      op,o2p,np,n2p,nop,ped,qsum,E1,E2,E3,u1,u2,u3,B1,B2,B3
    real(kind=rp),intent(out) :: nu_op,nu_o2p,nu_np,nu_n2p,nu_nop, &
      heat_i,heat_e,qei,qin,qen,kappa_i,kappa_e

    real(kind=rp),parameter :: ev2j = 1.602e-19_rp
    real(kind=rp) :: tr,ne,Eeq1,Eeq2,Eeq3

    tr = (tn + ti)/2
    ne = op + o2p + np + n2p + nop

    call collision_frequency(tr,o2,o1,he,n2, &
      nu_op,nu_o2p,nu_np,nu_n2p,nu_nop)
    call equivalent_efield(E1,E2,E3,u1,u2,u3,B1,B2,B3,Eeq1,Eeq2,Eeq3)
    heat_i = ion_heat(ped,Eeq1,Eeq2,Eeq3)
    heat_e = ev2j*qsum*electron_heat_coeff(ne/(o2 + o1 + n2))
    qei = ev2j*transfer_electron_ion(te,ne,op,o2p,np,n2p,nop)
    qin = ev2j*transfer_ion_neutral(tr,o2,o1,he,n2,op,o2p,np,n2p,nop)
    qen = ev2j*transfer_electron_neutral(tn,o2,o1,he,n2,ne,te)
    kappa_i = ev2j*ion_conductivity(op,ne)
    kappa_e = ev2j*electron_conductivity(te,ne,o2,o1,he,n2)

  endsubroutine calculate
!-----------------------------------------------------------------------
endmodule calculate_module