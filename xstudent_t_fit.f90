program simulate_student_t_fit
  use kind_mod, only: dp
  use random_mod, only: random_student_t
  use student_t_fit_mod, only: fit_student_t, fit_student_t_nu_fixed
  use basic_stats_mod, only: mean
  implicit none
  integer, parameter :: n = 10**5, nsim=10
  real(kind=dp), allocatable :: xdata(:), nu_try(:)
  real(kind=dp) :: mu0, s0, nu0
  real(kind=dp) :: mu_hat, s_hat, nu_hat
  integer :: isim, inu
  character (len=*), parameter :: fmt_cr = "(a20,':',*(1x,f12.6))", &
     fmt_cc = "(a20,':',*(1x,a12))"
  print*,"#obs:", n

  ! Set true parameters.
  mu0 = 0.5_dp    ! true location
  s0  = 1.2_dp    ! true scale
  nu0 = 3.0_dp    ! true degrees-of-freedom
  nu_try = real([1, 2, 3, 4, 5], kind=dp)

  allocate(xdata(n))
  print fmt_cc,"parameter", "mu", "s", "nu", "mean"
  print fmt_cr,"true", mu0, s0, nu0
  ! Simulate data: x = mu0 + s0 * t, where t ~ Student_t(nu0).
  ! Fit the Student t distribution to the simulated data.
  do isim=1,nsim
     xdata = mu0 + s0*random_student_t(n, nu0)
     call fit_student_t(xdata, mu_hat, s_hat, nu_hat)
     print fmt_cr,"estimated", mu_hat, s_hat, nu_hat, mean(xdata)
     do inu=1,size(nu_try)
        call fit_student_t_nu_fixed(xdata, nu_try(inu), mu_hat, s_hat)
        print fmt_cr,"estimated", mu_hat, s_hat, nu_try(inu), mean(xdata)
     end do
     print*
  end do
end program simulate_student_t_fit
