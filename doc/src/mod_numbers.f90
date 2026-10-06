module numbers 
  implicit none
  double precision, parameter :: ZERO=0.0d0, ONE=1.0d0, TWO=2.0d0, THREE=3.0d0, FOUR=4.0d0, FIVE=5.0d0
  double precision, parameter :: SIX=6.0d0, SEVEN=7.0d0, EIGHT=8.0d0, NINE=9.0d0, TEN=10.0d0
  double precision, parameter :: HALF=0.5d0, ONETHIRSD=1.0d0/3.0d0, TWOTHIRD=2.0d0/3.0d0, ONEANDHALF=1.5d0
  double precision, parameter :: QUARTER=0.25d0, THREEQUARTER=0.75d0
  double precision, parameter :: pi=4.d0 * datan(1.d0)
  double precision, parameter :: hr=1.d0*(3.d0**(-1.d0))
  double precision, parameter :: boh=5.d0*(3.d0**(-1.d0))
  double precision, parameter :: bs=1.d0*(6.d0**(-1.d0)) 
  double precision, parameter :: bih=2.d0*(3.d0**(-1.d0))

  contains
    function dfact(a)
      integer :: dfact
      integer :: i
      integer, intent(in) :: a
      dfact=1
      if (mod(a,2).eq.0) then !n is even
        do i=1,a/2
          dfact=dfact*(2*i)
        end do
      else 
        do i=1,(a+1)/2
          dfact=dfact*(2*i-1)
        end do
      end if
    end function  
end module numbers