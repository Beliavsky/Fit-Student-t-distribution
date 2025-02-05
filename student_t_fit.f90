module student_t_fit_mod
  use kind_mod, only: dp
  use nelder_mead_mod, only: nelder_mead
  use obj_fun_mod,   only: obj_fun, dist_type, xdata, nu0
  use constants_mod, only: pi
  implicit none
  private
  public :: fit_student_t, fit_student_t_nu_fixed

contains

  !-----------------------------------------------------------------
  ! Subroutine: fit_student_t
  !
  ! Fits a Student t distribution (location, scale, degrees-of-freedom)
  ! to the univariate data in the array DATA using maximum likelihood.
  ! The fitting is performed via the Nelder–Mead optimizer.
  !
  ! Input:
  !   data     - univariate data (array)
  ! nu_guess   - guess for degrees-of-freedom
  ! Output:
  !   mu_hat   - estimated location parameter
  !   s_hat    - estimated scale parameter
  !   nu_hat   - estimated degrees-of-freedom
  !-----------------------------------------------------------------
  subroutine fit_student_t(xx, mu_hat, s_hat, nu_hat, nu_guess)
    real(kind=dp), intent(in)  :: xx(:)
    real(kind=dp), intent(out) :: mu_hat, s_hat, nu_hat
    real(kind=dp), intent(in), optional :: nu_guess
    real(kind=dp) :: x0(3), xopt(3), fopt
    integer :: info, niter
    real(kind=dp), parameter :: tol = 1.0e-6, nu_guess_=5.0_dp
    integer, parameter :: max_iter = 1000
    ! Store the data in the module variable.
    xdata = xx
    ! Set the distribution type so that the global objective function
    ! (obj_fun in obj_fun_mod) will call neg_loglik_student_t.
    dist_type = "student_t"
    ! Use sample mean and standard deviation as initial guesses.
    x0(1) = sum(xx) / real(size(xx), dp)
    x0(2) = sqrt( sum((xx - x0(1))**2) / real(size(xx), dp) )
    if (present(nu_guess)) then
       x0(3) = nu_guess
    else
       x0(3) = nu_guess_
    end if
    call nelder_mead(x0, max_iter, tol, xopt, fopt, info, niter)
    mu_hat = xopt(1)
    s_hat  = xopt(2)
    nu_hat = xopt(3)
  end subroutine fit_student_t

  !-----------------------------------------------------------------
  ! Subroutine: fit_student_t_nu_fixed
  !
  ! Fits a Student t distribution (location, scale) with fixed nu
  ! to the univariate data in the array DATA using maximum likelihood.
  ! The fitting is performed via the Nelder–Mead optimizer.
  !
  ! Input:
  !   data     - univariate data (array)
  ! nu_fixed   - degrees-of-freedom
  ! Output:
  !   mu_hat   - estimated location parameter
  !   s_hat    - estimated scale parameter
  !   nu_hat   - estimated degrees-of-freedom
  !-----------------------------------------------------------------
  subroutine fit_student_t_nu_fixed(xx, nu_fixed, mu_hat, s_hat)
    real(kind=dp), intent(in)  :: xx(:)
    real(kind=dp), intent(in)  :: nu_fixed
    real(kind=dp), intent(out) :: mu_hat, s_hat
    real(kind=dp) :: x0(2), xopt(2), fopt
    integer :: info, niter
    real(kind=dp), parameter :: tol = 1.0e-6
    integer, parameter :: max_iter = 1000
    ! Store the data in the module variable.
    xdata = xx
    ! Set the distribution type so that the global objective function
    ! (obj_fun in obj_fun_mod) will call neg_loglik_student_t.
    dist_type = "student_t_nu_fixed"
    nu0 = nu_fixed
    ! Use sample mean and standard deviation as initial guesses.
    x0(1) = sum(xx) / real(size(xx), dp)
    x0(2) = sqrt( sum((xx - x0(1))**2) / real(size(xx), dp) )
    call nelder_mead(x0, max_iter, tol, xopt, fopt, info, niter)
    mu_hat = xopt(1)
    s_hat  = xopt(2)
  end subroutine fit_student_t_nu_fixed

end module student_t_fit_mod
