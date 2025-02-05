module obj_fun_mod
  use kind_mod, only: dp
  use constants_mod, only: pi, log_two, log_two_pi, sqrt_two
  use basic_stats_mod, only: variance
  implicit none
  private
  public :: obj_fun, xdata, dist_type, xdist, bad_nll, nu0

  ! Global data vector
  real(kind=dp), allocatable :: xdata(:)
  ! Global control variables:
  !   dist_type: selects the overall negative log likelihood to use.
  character(len=100) :: dist_type
  character(len=10)  :: xdist
  real(kind=dp)      :: nu0
  ! Penalty value returned if parameters are out-of-bounds.
  real(kind=dp), parameter :: bad_nll = 1.0e10_dp

contains

  !-----------------------------------------------------------------
  ! Function: obj_fun
  !
  ! This is the objective function to be minimized by the optimizer.
  ! It selects the appropriate negative log likelihood function based
  ! on the value of the global variable dist_type.
  !-----------------------------------------------------------------
  function obj_fun(x) result(y)
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: y
    select case (dist_type)
    case ("student_t")
       y = neg_loglik_student_t(x)
    case ("student_t_nu_fixed")
       y = neg_loglik_student_t_nu_fixed(x)
    case default
       y = bad_nll
    end select
  end function obj_fun

  !-----------------------------------------------------------------
  ! Function: neg_loglik_student_t
  !
  ! Computes the negative log likelihood for a Student t distribution with
  ! parameters:
  !    par(1) = mu      (location),
  !    par(2) = s       (scale; must be > 0),
  !    par(3) = nu      (degrees-of-freedom; must be > 0).
  !
  ! The density is assumed to be
  !
  !   f(x|mu,s,nu) = Gamma((nu+1)/2) / ( Gamma(nu/2)*sqrt(nu*pi)*s ) *
  !                  [1 + ((x-mu)/s)**2/nu]^{-(nu+1)/2}.
  !
  ! For any parameter value that violates s > 0 or nu > 0, a large penalty
  ! (bad_nll) is returned.
  !-----------------------------------------------------------------
  function neg_loglik_student_t(par) result(nll)
    real(kind=dp), intent(in) :: par(:)
    real(kind=dp) :: nll, mu, s, nu, tmp, term
    integer :: i, n

    n = size(xdata)
    mu = par(1)
    s  = par(2)
    nu = par(3)

    if (s <= 0.0_dp .or. nu <= 0.0_dp) then
       nll = bad_nll
       return
    end if

    nll = 0.0_dp
    do i = 1, n
       tmp = (xdata(i) - mu) / s
       term = 1.0_dp + (tmp**2)/nu
       nll = nll + ( - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) + 0.5_dp*log(nu*pi) + log(s)  &
                     + ((nu+1.0_dp)/2.0_dp)*log(term) )
    end do
  end function neg_loglik_student_t

  !-----------------------------------------------------------------
  ! Function: neg_loglik_student_t
  !
  ! Computes the negative log likelihood for a Student t distribution with
  ! parameters:
  !    par(1) = mu      (location),
  !    par(2) = s       (scale; must be > 0),
  !    par(3) = nu      (degrees-of-freedom; must be > 0).
  !
  ! The density is assumed to be
  !
  !   f(x|mu,s,nu) = Gamma((nu+1)/2) / ( Gamma(nu/2)*sqrt(nu*pi)*s ) *
  !                  [1 + ((x-mu)/s)**2/nu]^{-(nu+1)/2}.
  !
  ! For any parameter value that violates s > 0 or nu > 0, a large penalty
  ! (bad_nll) is returned.
  !-----------------------------------------------------------------
  function neg_loglik_student_t_nu_fixed(par) result(nll)
    real(kind=dp), intent(in) :: par(:)
    real(kind=dp) :: nll, mu, s, nu, tmp, term
    integer :: i, n

    n = size(xdata)
    mu = par(1)
    s  = par(2)
    nu = nu0
    if (s <= 0.0_dp .or. nu <= 0.0_dp) then
       nll = bad_nll
       return
    end if

    nll = 0.0_dp
    do i = 1, n
       tmp = (xdata(i) - mu) / s
       term = 1.0_dp + (tmp**2)/nu
       nll = nll + ( - log_gamma((nu+1.0_dp)/2.0_dp) + log_gamma(nu/2.0_dp) + 0.5_dp*log(nu*pi) + log(s)  &
                     + ((nu+1.0_dp)/2.0_dp)*log(term) )
    end do
  end function neg_loglik_student_t_nu_fixed

end module obj_fun_mod
