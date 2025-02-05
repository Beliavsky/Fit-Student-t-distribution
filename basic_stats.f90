module basic_stats_mod
use kind_mod, only: dp
implicit none
private
public :: mean, variance
contains

  !------------------------------------------------------------
  ! Function to compute the sample variance of an array.
  !------------------------------------------------------------
  function variance(x) result(var)
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: var, m
    integer :: n
    n = size(x)
    m = sum(x) / n
    var = sum((x - m)**2) / (n-1)
  end function variance

function mean(x) result(xmean)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: xmean
xmean = sum(x)/max(1,size(x))
end function mean

end module basic_stats_mod