module types

  use :: prcn_mod
  implicit none

  type,public :: cart_position
    real(wp) :: x,y,z
  end type cart_position

  type,public :: decomp
    logical :: Success
  end type decomp

  type :: Ptr2Real
    real(wp),pointer :: p
  end type Ptr2Real

end module types
