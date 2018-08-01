module cmn_tp_mod

  use :: prcn_mod

  ! Components of displacement between two beads
  type :: dis
    real(wp) :: x,y,z
    real(wp) :: mag,mag2
    real(wp) :: riy,rjy,yim
    real(wp) :: magim,mag2im
  end type

end module cmn_tp_mod
