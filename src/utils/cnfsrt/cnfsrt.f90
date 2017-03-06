
  ! This routine is written by Vydia Venkataramani and modified by Amir Saadat.
  ! It is for sorting the ensemble to different configurational classes:
  ! Description of diffrent configuration types:
  ! Type 1   :   Folds(2/3 or 2/5) 
  ! Type 2   :   Half-Dumbells 
  ! Type 3   :   Kinks
  ! Type 4   :   Dumbells 
  ! Type 5   :   Coils
  ! Type 6   :   Extended
  subroutine conf_sort(q,rvmrc,nseg,nbead,cnf_tp)

    use :: inp_dlt, only: npchain,qmax,residx
    use :: arry_mod, only: print_vector,print_matrix

    integer,intent(inout) :: nseg,nbead,cnf_tp(:)
    real(wp),intent(in),target :: q(:,:),rvmrc(:,:)
    real(wp),pointer :: qPx(:),qPy(:),qPz(:)
    real(wp),pointer :: RPx(:),RPy(:),RPz(:)
!    real(wp) :: pconfig(6)
!    integer :: iconfig(6)
    real(wp) :: rmax,xi,yi,zi,xj,yj,zj,dist,xmin,ymin,zmin
    real(wp) :: xmax,ymax,zmax,xr,yr,zr,xcap,ycap,zcap,rmin
    real(wp) :: xproj,yproj,zproj,xpres,xnext,dia,radius
    real(wp) :: bead_left,bead_right,avg,avgsum
    integer :: i,ntype,maxlength,ichain,ibead,jbead,ilen
    integer :: min,max,maxbright,maxbright1,low,last,info,info1
    integer :: k,ij,j,ibcount,sumbright,lcount,left,markerc,ik
    integer :: lend,marker1,marker2,markerl,markerr,icount

    if (.not.allocated(x)) allocate(x(nbead))
    if (.not.allocated(y)) allocate(y(nbead))
    if (.not.allocated(z)) allocate(z(nbead))
    if (.not.allocated(r)) allocate(r(nbead))
    maxlength=int(qmax*nseg*residx)+1
    if (.not.allocated(ib))  allocate(ib(maxlength))
    if (.not.allocated(ib1)) allocate(ib1(maxlength))
    if (.not.allocated(icount3)) allocate(icount3(maxlength))

    do i=1, 6
!      iconfig(i)=0
!      pconfig(i)=0._wp
    end do

    do ichain=1, npchain
      qPx => q(1:3*nseg-2:3,ichain)
      qPy => q(2:3*nseg-1:3,ichain) 
      qPz => q(3:3*nseg:3,ichain)
      RPx => rvmrc(1:3*nbead-2:3,ichain)
      RPy => rvmrc(2:3*nbead-1:3,ichain) 
      RPz => rvmrc(3:3*nbead:3,ichain)
      ntype=0
      ! Calculating the position of each bead
!      x(1)=0._wp
!      y(1)=0._wp
!      z(1)=0._wp
!      do ibead=2, nbead
!        x(ibead)=x(ibead-1)+qPx(ibead-1)
!        y(ibead)=y(ibead-1)+qPy(ibead-1)
!        z(ibead)=z(ibead-1)+qPz(ibead-1)
!      end do          

      ib=0
      ib1=0
!      ! Calculating the largest distance between the beads
!      rmax=0._wp
!      do ibead=1, nbead-1
!        xi=x(ibead)
!        yi=y(ibead)
!        zi=z(ibead)
!        do jbead=1, nbead
!          xj=x(jbead)
!          yj=y(jbead)
!          zj=z(jbead)
!          dist=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
!
!          if (dist > rmax) then
!            min=ibead
!            max=jbead
!            rmax=dist
!          end if
!        end do ! jbead
!      end do ! ibead
!
!      xmin=x(min)
!      ymin=y(min)
!      zmin=z(min)
!
!      xmax=x(max)
!      ymax=y(max)
!      zmax=z(max)
!
!      ! Moving beads w.r.t min
!      do ibead=1, nbead
!         x(ibead)=x(ibead)-xmin
!         y(ibead)=y(ibead)-ymin
!         z(ibead)=z(ibead)-zmin
!      end do

      ! Calculating the largest distance between the beads
      rmax=0._wp
      do ibead=1, nbead-1
        xi=RPx(ibead)
        yi=RPy(ibead)
        zi=RPz(ibead)
        do jbead=1, nbead
          xj=RPx(jbead)
          yj=RPy(jbead)
          zj=RPz(jbead)
          dist=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)

          if (dist > rmax) then
            min=ibead
            max=jbead
            rmax=dist
          end if
        end do ! jbead
      end do ! ibead

      xmin=RPx(min)
      ymin=RPy(min)
      zmin=RPz(min)

      xmax=RPx(max)
      ymax=RPy(max)
      zmax=RPz(max)

      do ibead=1, nbead
        x(ibead)=RPx(ibead)-xmin
        y(ibead)=RPy(ibead)-ymin
        z(ibead)=RPz(ibead)-zmin
      end do
!call print_vector(x,'x2')
!call print_vector(y,'y2')
!call print_vector(z,'z2')

      ! Calculating the unit vector along the longest length
      xr=xmax-xmin
      yr=ymax-ymin
      zr=zmax-zmin

      xcap=xr/rmax
      ycap=yr/rmax
      zcap=zr/rmax

      ! Finding the projection of each bead onto the longest length
      rmin=0._wp
      do ibead=1, nbead
        r(ibead)=0._wp
        xproj=x(ibead)*xcap
        yproj=y(ibead)*ycap
        zproj=z(ibead)*zcap
        r(ibead)=(xproj+yproj+zproj)*residx
        if (r(ibead) < rmin) then
          rmin=r(ibead)
        end if
      end do

      ! Repositioning the beads such that the minimum r = 0
      do ibead=1, nbead
        r(ibead)=r(ibead)-rmin
      end do

      ! Accounting for bond overlaps
      maxbright=int(rmax*residx)+1
!print *,'rmax',rmax
!print *,'maxbright',maxbright
      if (r(1) > r(2)) then     
        low=int(r(2))+1
        last=int(r(1))+1
        info=0
      else
        low=int(r(1))+1 
        last=int(r(2))+1
        info=1
      end if

      do ilen=low, last
        ib(ilen)=ib(ilen)+1
      end do

      do ibead=2, nbead-1
        xpres=r(ibead)
        xnext=r(ibead+1)
        if (xpres < xnext) then
          info1 = 0
        else
          info1 = 1
        end if

        if ((info == 0).and.(info1 == 0)) then
          low=int(xpres)+2
          last=int(xnext)+1
        end if

        if ((info == 0).and.(info1 == 1)) then
          low=int(xnext)+1
          last=int(xpres)+1
        end if

        if ((info == 1).and.(info1 == 1)) then
          low=int(xnext)+1
          last=int(xpres)
        end if

        if ((info == 1).and.(info1 == 0)) then
          low=int(xpres)+1
          last=int(xnext)+1
        end if

        do ilen=low, last
          ib(ilen)=ib(ilen)+1
        end do
        info=info1
      end do

      ! Accounting for bead overlap

      dia=1._wp*residx ! parameter set based on your assumption
!      dia=qmax*residx ! parameter set based on your assumption
      radius=dia/2
           
      do i=1, maxbright
        icount3(i)=0._wp
      end do
      
      do ibead=1, nbead
        bead_left=r(ibead)-radius
        bead_right=r(ibead)+radius
        if (bead_left < 0._wp) then
          bead_left = 0._wp
        end if

        if (bead_right > real(maxbright,kind=wp)) then
          bead_right=real(maxbright,kind=wp)
        end if
        low=int(bead_left)+1
        last=int(bead_right)+1
        do ilen=low, last
          ib(ilen)=ib(ilen)+1
          icount3(ilen)=1
        end do
      end do

      do ilen=1, maxbright
        if (icount3(ilen) == 1) then
          ib(ilen)=ib(ilen)-1
        end if
      end do

      ! Averaging over configurations
!call print_vector(ib(1:maxbright),'ib1')

      if (mod(maxbright,residx) == 0) then
        maxbright1=maxbright/residx
      else
        maxbright1=int(maxbright/residx)+1
      end if

      k=0
      ij=1
      do i=1, maxbright1
        j=1
        ibcount=0
        do while ((j <= residx).and.(ij <= maxbright))
          ibcount=ibcount+ib(ij)
          j=j+1
          ij=ij+1
        end do
        if (j == residx) then
          avg=real(ibcount,kind=wp)/residx
        else
          avg=real(ibcount,kind=wp)/(j-1)
        end if

        if ((real(int(avg),kind=wp)+0.5_wp) < avg) THEN
          ib1(i)=int(avg)+1
        else
          ib1(i)=int(avg)
        end if
        k=k+1
      end do

      do i=1, maxbright
        ib(i)=0
      end do
      maxbright=maxbright1
      do i=1, maxbright
        ib(i)=ib1(i)
      end do

!call print_vector(ib(1:maxbright),'ib2')
      ! Checking for the minimum resolution

!      if (rmax <= 15._wp) then
      if (rmax <= 0.05_wp*qmax*nseg) then
!        iconfig(5)=iconfig(5)+1
        ntype=5
      else

        ! **Assigning for 2/3, 2/5 folds, coil and half dumbell**

        if (((ib(1) == 1).and.(ib(maxbright) > 1)).or. &
            ((ib(maxbright) == 1).and.(ib(1) > 1))) then
           
          if ((ib(1) == 1).and.(ib(maxbright) > 1)) then
             j=1
             i=maxbright
          end if
          if ((ib(maxbright) == 1).and.(ib(1) > 1)) then
            j=0
            i=1
          end if
           
          icount = 0   
          sumbright = 0
          do while (ib(i) > 1)
            icount=icount+1
            sumbright=sumbright+ib(i)
            if (j == 1) then
              i=i-1
            else
              i=i+1
            end if  
          end do
          avgsum=sumbright/(1._wp*icount)
          lcount=0
          if (j == 1) then
            k=1
          else
            k=maxbright
          end if
          do while (ib(k) == 1)
            lcount = lcount + 1
            if (j == 1) then
              k=k+1
            else
              k=k-1
            end if
          end do
          if (j == 1) then
            left=k
            last=i
          else
            left=i
            last=k 
          end if
          
          markerc=0
          do ik=left, last
            if (ib(ik) > (int(avgsum))) then
              markerc=1
              lend=ik
            end if
          end do

          if (markerc == 1) then
!            iconfig(3)=iconfig(3)+1
            ntype=3
          else
            icount=maxbright-lcount
            if ((0.25_wp <= ((1._wp*icount)/maxbright))) then
!              iconfig(1)=iconfig(1)+1
              ntype=1
            end if
            if (((1._wp*icount)/maxbright) < 0.25_wp) then
!              iconfig(2)=iconfig(2)+1
              ntype=2
            end if
          end if
        end if

        ! **Assigning for kink or fully extended**
              
        if ((ib(1) == 1).and.(ib(maxbright) == 1)) THEN
          marker1=0
          do i=2, (maxbright-1)
            if (ib(i) == 1) then
              marker1=marker1+1
            end if
          end do
          if (marker1 == (maxbright-2)) then
!            iconfig(6)=iconfig(6)+1
            ntype=6
          else
!            iconfig(3)=iconfig(3)+1
            ntype=3
          end if
        end if

        ! **Assigning for coil or dumbell**
              
        if ((ib(1) > 1).and.(ib(maxbright) > 1)) THEN
          marker2=0
          markerl=2
          do while ((ib(markerl) > 1).and.(markerl < maxbright))
            markerl=markerl+1
          end do
          markerr=maxbright-1
          do while ((ib(markerr) > 1).and.(markerr > 1))
            markerr=markerr-1
          end do
          if ((markerl == (maxbright)).and.(markerr == 1)) then
!            iconfig(5)=iconfig(5)+1
            ntype=5
          else
            do i=markerl, markerr
              if (ib(i) > 1) then
                marker2=1
              end if
            end do
            if (marker2 /= 1) then
!              iconfig(4)=iconfig(4)+1
              ntype=4
            else
!              iconfig(5)=iconfig(5)+1
              ntype=5
            end if
          end if           
        end if
              
      end if

      cnf_tp(ichain)=ntype
                
    end do ! ichain

!    do i=1, 6
!      pconfig(i)=real(iconfig(i),kind=wp)/nchain
!    end do

  end subroutine conf_sort
