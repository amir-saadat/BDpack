  subroutine update_bendforce(Fphi,nchain,nbead,nchian_cmb,nseg_cmb,&
                             nseg_cmbbb,nseg_cmbar,add_cmb,ForceLaw,Qt,
  
  
    use :: arry_mod, only: print_vector
    !use :: conv_mod, only: Bbar_vals,Bbar_cols,Bbar_rowInd
    use :: flow_mod, only: FlowType
    use :: force_smdlt, only: Fphi
    
    class(sprforce),intent(inout) :: this
    real(wp),intent(in) :: Rbx(:)
    real(wp),intent(in) :: Rby(:)
    real(wp),intent(in) :: Rbz(:)
    real(wp),intent(in) :: bs(3),invbs(3)
    real(wp),intent(out) :: Fbnd(:)
    integer,intent(in) :: itime,ntotseg,ntotsegx3,ntotbead,ntotbeadx3
	integer,intent(in) :: nchain,nseg,nbead,nseg_cmbar !not needed here 
    real(wp),intent(in) :: Qt(:)
    integer :: its,ich,osb,oss,is
    real(wp) :: qx,qy,qz,qsq,Ftmp,qytmp

    integer,intent(in) :: itime
    real(wp) :: thta(-1:1),cost(-1:1),thtal,thtar,costl,costr
    real(wp) :: thta_s,cost_s
    real(wp) :: qtmp(3,-2:1),qmg(-2:1),ehat(3,-2:1)
    real(wp) :: qtmpl(3),qtmpr(3),qmgl,qmgr,ehatl(3),ehatr(3)
    integer  :: nbead,ibead,os,nbead_cmb,osl,iarm

    allocate(Fbnd(ntotbeadx3))	
    Fbnd=0	
	if (nbead=>3 .and. nchian =/0) then
     do ichian=1, nchian
	 
       Osb1=(ichian-1)*nbead
	   OsS1=(ichian-1)*nSeg
	   
       do ibead=3, nbead
	   
	      osS=OsS1+(ibead-1)
	      osb=Osb1+ ibead
		  
	      qtmp(:,-1)=Qt((osS-1)*3+1:(osS-1)*3+3)
	      qtmp(:,-2)=Qt((osS-2)*3+1:(osS-2)*3+3)
	      qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2))) 
	      qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) 
	      ehat(:,-2)=qtmp(:,-2)/qmg(-2)
	      ehat(:,-1)=qtmp(:,-1)/qmg(-1) 
	      thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
	      cost(-1)=cos(thta(-1))
	      !!Fi
	      Fbnd(osb*3-2:osb*3+0)=Fbnd(osb*3-2:osb*3)+&
		               WLC_C*(1/qmg(-1))*(ehat(:,-2)-cost(-1)*ehat(:,-1))
	      !!Fi-1
	      Fbnd((osb-1)*3-2:(osb-1)*3+0)=Fbnd((osb-1)*3-2:(osb-1)*3)+&
		              WLC_C*( ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))-&
                              ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2)) )
	      !!Fi-2
	      Fbnd((osb-2)*3-2:(osb-2)*3+0)=Fbnd((osb-2)*3-2:(osb-2)*3)+&
		              WLC_C*(1/qmg(-2))*( cost(-1)*ehat(:,-2)-ehat(:,-1))
       end do !ibead
     end do !ichain
	end if


    if (add_cmb) then
	
     do ichian=1, nchian_cmb
	 
      Osb1=nchian*nbead+(ichian-1)*nbead_cmb
	  OsS1=nchian*nSeg+(ichian-1)*nSeg_cmb
	  ! Loop over backbone
       do ibead=3, nbead_cmbbb
	   
	      osS=OsS1+(ibead-1)
	      osb=Osb1+ ibead
		  
	      qtmp(:,-1)=Qt((osS-1)*3+1:(osS-1)*3+3)
	      qtmp(:,-2)=Qt((osS-2)*3+1:(osS-2)*3+3)
	      qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2))) 
	      qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) 
	      ehat(:,-2)=qtmp(:,-2)/qmg(-2)
	      ehat(:,-1)=qtmp(:,-1)/qmg(-1) 
	      thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
	      cost(-1)=cos(thta(-1))
	      !!Fi
	      Fbnd(osb*3-2:osb*3+0)=Fbnd(osb*3-2:osb*3)+&
		               WLC_C*(1/qmg(-1))*(ehat(:,-2)-cost(-1)*ehat(:,-1))
	      !!Fi-1
	      Fbnd((osb-1)*3-2:(osb-1)*3+0)=Fbnd((osb-1)*3-2:(osb-1)*3)+&
		              WLC_C*( ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))-&
                              ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2)) )
	      !!Fi-2
	      Fbnd((osb-2)*3-2:(osb-2)*3+0)=Fbnd((osb-2)*3-2:(osb-2)*3)+&
		              WLC_C*(1/qmg(-2))*( cost(-1)*ehat(:,-2)-ehat(:,-1))

       end do !ibead_cmbbb
	 
	   ! Loop over Arms
	   ! Bead[end of BB]
	   Os1 = nchian*nbead + (ichian-1)*nbead_cmb + nbead_cmbbb
	   !Segment Qt[end of BB]
	   OsS1= nchian*nSeg + (ichian-1)*nSeg_cmb  + nSeg_cmbbb
      
if     (seg_cmbar==1) then
      do iarm=1, Na
!Backbone _Arm
	    !Segment of the backbone connected to the arm
		!Ia: Backbone bead place of arm
        osSbb=nchian*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)  
		!bead id for iarm start  !Ia: Backbone bead place of arm
		oslbbb=nchian*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))  
!Seg-Arm
		osS=OsS1+(iarm-1)*nseg_cmbar  
        ! osb+1 =Bead# start arm iarm
		osb=os1+(iarm-1)*(seg_cmbar)   ! nbead_arm = seg_cmbarm
!-------------------------------------
! Effect of first bead of arm on the Backbone.
! Left and right segments 
! ...O-L-O-R-O...
!        | arm Seg 1
!        O arm bead 1
! -------------------------------
		qtmpl(:)= q(osSbb*3-2:osSbb*3)
		qtmpr(:)=-q(osSbb*3+1:osSbb*3+3)
		qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
		qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
		ehatl(:)=qtmpl(:)/qmgl
		ehatr(:)=qtmpr(:)/qmgr
		
		qtmp(:,1)=q(osS*3+1:osS*3+3)  ! 1st segment of the arm
		qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1))) !1st_seg arm
		ehat(:,1)=qtmp(:,1)/qmg(1)  !1st_seg arm
		thtal=acos(dot(qtmpl,qtmp(:,1))/(qmgl*qmg(1)))
		thtar=acos(dot(qtmpr,qtmp(:,1))/(qmgr*qmg(1)))
		costl=cos(thtal)
		costr=cos(thtar)
		! force on BB Bead Left  n-1
		Fbnd(oslbbb*3-5:oslbbb*3-3)=Fbnd(oslbbb*3-2:oslbbb*3)+&
                                        WLC_C/qmgl*(costl*ehatl(:)-ehat(:,1))
		! force on BB Bead Right n+1
		Fbnd(oslbbb*3+1:oslbbb*3+3)=Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                        WLC_C/qmgr*(costr*ehatr(:)-ehat(:,1))
        ! force on BB Bead Right n
		Fbnd(oslbbb*3-2:oslbbb*3+0)=Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C*&
              (  (ehat(:,1)*(1/qmgl  +costl/qmg(1)) -ehatl(:)*(1/qmg(1)+costl/qmgl))+&
                 (ehat(:,1)*(1/qmgr  +costr/qmg(1)) -ehatr(:)*(1/qmg(1)+costr/qmgr)) )
		! 1st bead of the arm 2@F(v,v-1)
		Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( 1/qmg(1)*(ehatr(:)-costr*ehat(:,1))+&
		                                  1/qmg(1)*(ehatl(:)-costl*ehat(:,1)) )
		end do
		
elseif (seg_cmbar==2) then
      do iarm=1, Na
!Backbone _Arm
	    !Segment of the backbone connected to the arm
		!Ia: Backbone bead place of arm
        osSbb=nchian*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)  
		!bead id for iarm start  !Ia: Backbone bead place of arm
		oslbbb=nchian*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))  
!Seg-Arm
		osS=OsS1+(iarm-1)*seg_cmbar  
        ! osb/osS +1 = First Beadtot#/Segtot# start arm iarm   
		osb=os1+(iarm-1)*(seg_cmbar)   ! nbead_arm = seg_cmbarm
!-------------------------------------
! Effect of first bead of arm on the Backbone.
! Left and right segments 
! ...O-L-O-R-O...
!        | arm Seg 1
!        O arm bead 1
! -------------------------------
		qtmpl(:)= q(osSbb*3-2:osSbb*3)
		qtmpr(:)=-q(osSbb*3+1:osSbb*3+3)
		qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
		qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
		ehatl(:)=qtmpl(:)/qmgl
		ehatr(:)=qtmpr(:)/qmgr
		
		qtmp(:,1)=q(osS*3+1:osS*3+3)  ! First segment of the arm
		qtmp(:,2)=q((osS+1)*3+1:(osS+1)*3+3)
		qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1))) !1st_seg arm
		qmg(2)=sqrt(dot(qtmp(:,2),qtmp(:,2)))
		ehat(:,1)=qtmp(:,1)/qmg(1)  !1st_seg arm
		ehat(:,2)=qtmp(:,2)/qmg(2)
		thta(1)=acos(dot(qtmp(:,1),qtmp(:,2))/(qmg(1)*qmg(2))) !Angle Seg1&Seg2
		cost(1)=cos(thta(1))
		! force on BB Bead Left  n-1
		Fbnd(oslbbb*3-5:oslbbb*3-3)=Fbnd(oslbbb*3-2:oslbbb*3)+&
                                        WLC_C/qmgl*(costl*ehatl(:)-ehat(:,1))
            ! force on BB Bead Right n+1
		Fbnd(oslbbb*3+1:oslbbb*3+3)=Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                         WLC_C/qmgr*(costr*ehatr(:)-ehat(:,1))
			! force on BB Bead Right n
		Fbnd(oslbbb*3-2:oslbbb*3+0)=Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C* &
                ( (ehat(:,1)*(1/qmgl  +costl/qmg(1)) -ehatl(:)*(1/qmg(1)+costl/qmgl))+&
                  (ehat(:,1)*(1/qmgr  +costr/qmg(1)) -ehatr(:)*(1/qmg(1)+costr/qmgr))+&
                  (1/qmg(1))*(cost(2)*ehat(:,1)-ehat(:,2)) )

		! 1st bead of the arm 2@F(v,v-1)+F(v,v)
		Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( 1/qmg(1)*(ehatr(:)-costr*ehat(:,1))+&
		                                  1/qmg(1)*(ehatl(:)-costl*ehat(:,1))+& 
										  (ehat(:, 2)*(1/qmg(1)+cost(1)/qmg(2))-&
                                          ehat(:, 1)*(1/qmg(2)+cost(1)/qmg(1))) )
		! 2ed Bead F(v,v-1)
		Fbnd((osb+1)*3+1:(osb+1)*3+3)=WLC_C*(0+& 
										  1/qmg(2)*(ehat(:,1)-cost(1)*ehat(:,1))

      end do


elseif (seg_cmbar >= 3)) then
      do iarm=1, Na
!Backbone _Arm
	    !Segment of the backbone connected to the arm
		!Ia: Backbone bead place of arm
        osSbb=nchian*nSeg+(ichain-1)*nSeg_cmb+(Ia(iarm+1)-1)  
		!bead id for iarm start  !Ia: Backbone bead place of arm
		oslbbb=nchian*nbead +(ichain-1)*(nSeg_cmb+1)+(Ia(iarm+1))  
!Seg-Arm
		osS=OsS1+(iarm-1)*seg_cmbar    ! +1 =1st Seg/bead iarm
		osb=os1+(iarm-1)*(seg_cmbar)   ! nbead_arm = seg_cmbarm
!-------------------------------------
! Effect of first bead of arm on the Backbone.
! Left and right segments 
! ...O-L-O-R-O...
!        | arm Seg 1
!        O arm bead 1
! -------------------------------
		qtmpl(:)= q(osSbb*3-2:osSbb*3)
		qtmpr(:)=-q(osSbb*3+1:osSbb*3+3)
		qmgl=sqrt(dot(qtmpl(:),qtmpl(:)))
		qmgr=sqrt(dot(qtmpr(:),qtmpr(:)))
		ehatl(:)=qtmpl(:)/qmgl
		ehatr(:)=qtmpr(:)/qmgr
		
        qtmp(:,1)=q(osS*3+1:osS*3+3)  ! 1st segment of the arm
		qtmp(:,2)=q((osS+1)*3+1:(osS+1)*3+3) !! 2ed segment of the arm
		
		qmg(1)=sqrt(dot(qtmp(:,1),qtmp(:,1))) !1st_seg arm
		ehat(:,1)=qtmp(:,1)/qmg(1)  !1st_seg arm
		qmg(2)=sqrt(dot(qtmp(:,2),qtmp(:,2)))
		ehat(:,2)=qtmp(:,2)/qmg(2)
		
		thtal=acos(dot(qtmpl,qtmp(:,1))/(qmgl*qmg(1)))
		thtar=acos(dot(qtmpr,qtmp(:,1))/(qmgr*qmg(1)))
		costl=cos(thtal)
		costr=cos(thtar)
		
		thta(1)=acos(dot(qtmp(:,1),qtmp(:,2))/(qmg(1)*qmg(2)))
        cost(1)=cos(thta(1))
		
        ! force on BB Bead Left  n-1
        Fbnd(oslbbb*3-5:oslbbb*3-3)=Fbnd(oslbbb*3-2:oslbbb*3)+&
                                        WLC_C/qmgl*(costl*ehatl(:)-ehat(:,1))
        ! force on BB Bead Right n+1
        Fbnd(oslbbb*3+1:oslbbb*3+3)=Fbnd(oslbbb*3+1:oslbbb*3+3)+&
                                         WLC_C/qmgr*(costr*ehatr(:)-ehat(:,1))
        ! force on BB Bead Right n
        Fbnd(oslbbb*3-2:oslbbb*3+0)=Fbnd(oslbbb*3-2:oslbbb*3)+WLC_C* &
                ( (ehat(:,1)*(1/qmgl  +costl/qmg(1)) -ehatl(:)*(1/qmg(1)+costl/qmgl))+&
                  (ehat(:,1)*(1/qmgr  +costr/qmg(1)) -ehatr(:)*(1/qmg(1)+costr/qmgr))+&
                  (1/qmg(1))*(cost(2)*ehat(:,1)-ehat(:,2)) )
        !1st bead of the arm 2@F(v,v-1)+F(v,v) !+  F(v,v+1) on the loop
		Fbnd((osb)*3+1:(osb)*3+3)=WLC_C*( 1/qmg(1)*(ehatr(:)-costr*ehat(:,1))+&
		                                  1/qmg(1)*(ehatl(:)-costl*ehat(:,1))+& 
										  (ehat(:, 2)*(1/qmg(1)+cost(1)/qmg(2))-&
                                           ehat(:, 1)*(1/qmg(2)+cost(1)/qmg(1))) )
		! 2ed Bead F(v,v-1) !+F(v,v)+ F(v,v+1) on the loop
		Fbnd((osb+1)*3+1:(osb+1)*3+3)=WLC_C*(1/qmg(2)*(ehat(:,1)-cost(1)*ehat(:,2))
		
		
        do ibead_arm=3,seg_cmbar
		   osS=OsS1+(iarm-1)*seg_cmbar +ibead_arm
		   osb=os1+(iarm-1)*(seg_cmbar)+ibead_arm ! seg_cmbarm=nbead_arm
		   qtmp(:,0)=q((osS-1)*3+1:(osS-1)*3+3)
		   qtmp(:,-1)=q((osS-2)*3+1:(osS-2)*3+3)
		   qtmp(:,-2)=q((osS-3)*3+1:(osS-3)*3+3)
		   qmg(-2)=sqrt(dot(qtmp(:,-2),qtmp(:,-2))) !qmg(-1)
           qmg(-1)=sqrt(dot(qtmp(:,-1),qtmp(:,-1))) !qmg( 0)
           qmg( 0)=sqrt(dot(qtmp(:,0),qtmp(:,0)))
		   ehat(:,-2)=qtmp(:,-2)/qmg(-2)
		   ehat(:,-1)=qtmp(:,-1)/qmg(-1) 
           ehat(:,0)=qtmp(:,0)/qmg(0)
           !ehat(:,-2)=ehat(:,-1)
           !ehat(:,-1)=ehat(:, 0)
		   thta(-1)=acos(dot(qtmp(:,-1),qtmp(:,-2))/(qmg(-1)*qmg(-2)))
           cost(-1)=cos(thta(-1))

		   
		   !Fi
		   Fbnd(osb*3-2:osb*3+0)=Fbnd(osb*3-2:osb*3)+&
		               WLC_C*(1/qmg(-1))*(ehat(:,-2)-cost(-1)*ehat(:,-1))
		   !Fi-1
		   Fbnd((osb-1)*3-2:(osb-1)*3+0)=Fbnd((osb-1)*3-2:(osb-1)*3)+&
		              WLC_C*( ehat(:,-1)*(1/qmg(-2)+cost(-1)/qmg(-1))-&
                              ehat(:,-2)*(1/qmg(-1)+cost(-1)/qmg(-2)) )
		   !Fi-2
		   Fbnd((osb-2)*3-2:(osb-2)*3+0)=Fbnd((osb-2)*3-2:(osb-2)*3)+&
		              WLC_C*(1/qmg(-2))*( cost(-1)*ehat(:,-2)-ehat(:,-1))
	
		end do
      end do

	end do !ichain_comb
    end if !add comb

end subroutine update_bendforce