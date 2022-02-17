module Precision
  implicit none  
  integer, parameter :: double = kind(1.0d0)
  real(kind=double)  :: zero=0.0d0
  real(kind=double)  :: one=1.0d0
  public :: double, zero
end module Precision
