!**** fgvolumecalc.f90 ****

module volume_calc_mod

	implicit none
	contains

    ! bottom, top, left, right are in units of drfine and dzfine.
    ! points in the geometry are in units of drfine and dzfine

    real(8) function corrected_cell_volume( bottom, top, left, &
                        right, drfine, dzfine, polygons, npoly, volcorrected, cellcounter, flag)
		use types_mod
        implicit none

        integer npoly
        type(polygon_type) polygons(npoly)
        real(8) drfine, dzfine
        real(8) :: pi = 3.14159265358979
        real(8) bottom, top, left, right, tmp
	logical volcorrected
	integer cellcounter, flag

!  array of volumes computed every time a piece of cell is sliced off
!  by a segment. The size of the array is rather arbitrary, but must be
!  at least 2. More segment crossings than that is a problem.
        real(8) volume(4),volume0,volsum

	! vector variables:
	real(8) wx, wz, gx, gz
		
	! used to pass in the point on the segment into "volcell"
	real(8) xseg, zseg, dxstep

	! temporary variables for calculation of vectors
	real(8) r1, r2, z1, z2
	real(8) den

	integer i, j, k, crossings

        ! If inter_segment=.true. then the segment being worked on intersects
        ! the cell
	logical inter_segment

        ! If inter_total=.false. then no segments intersected this cell;
        ! in this case just use the standard rectangular volume formula
        ! if it is .true. then an intersection with at least one cell took
        ! place; use the shaved-off volume algorithm. This works by computing
        ! the volume inside the computing region for each intersection segment
        ! and adding them up to get the total cell volume in the computing
        ! region, i.e., in the region that is not inside metal.
	logical inter_total


        type(point_type) hitsbottom, hitstop, hitsleft, hitsright

        ! This routine computes the volumes of cells, including the correction
        ! required for those cell that have had corners or even big chunks shaved
        ! off by polygon segments that intersect them. Cell information is passed in
        ! and the cell is checked against every polygon segment to see if an
        ! intersection takes place. If it does the corrected volume is computed;
        ! if more than one intersection takes place the volumes are accumulated;
        ! if no intersection takes place corrected_cell_volume is loaded with the
        ! full rectangular volume.

        ! initialize inter_total
	inter_total = .false.

        ! start the corrected volume out at 0. Everytime a segment cuts through the cell and
        ! we calculate the volume on the outside of the cut and add it to volume. This allows
        ! for a proper calculation when two segments cut through a cell.
        volume0= pi*(right*dzfine - left*dzfine)*((top*drfine)**2 - (bottom*drfine)**2)
        do i=1,4
           volume(i)=0.d0
        end do

! set the number of segment crossings to zero for this cell at the start
        crossings=0

        ! loop over the polygons
        do i = 1,npoly

            !loop over the segments: j is point 1 and k is point 2.
            do j = 1,polygons(i)%numpoints

                ! Define k as the next point. If we have reached the last point
                ! use the first point in the polygon as the second point.
                if( j == polygons(i)%numpoints ) then
                    k = 1
                else
                    k = j + 1
                end if


                ! only process the last segment if this is a closed polygon
                if( j < polygons(i)%numpoints .or. polygons(i)%polytype == 1 ) then

                    ! Set inter_segment to .true. if this segment intersects the cell
                    ! and to .false. if it doesn't. Note that this logical function
                    ! also returns the important quantities hitsbottom, hitstop, hitsleft
                    ! and hitsleft which are necessary to compute the volume of a shaved
                    ! off region.
                    inter_segment = doesintersect( polygons(i)%points(j), &
                                                polygons(i)%points(k), &
                                                bottom, top, left, right, &
                                                hitsbottom, hitstop, hitsleft, hitsright, cellcounter, flag)

                    if( inter_segment == .true. ) then
                        ! calculate the shaved-off volume with the function volcell
                        ! after the quantities this routine needs are calculated

			! wx, wz -- unit vectors tangent to the segment
			! gx, gz -- unit vectors pointing INTO the metal
			! xseg, zseg -- a point on the segment (inside the collision cell)
			! dxstep -- defines tolerance for real equality check (say...drfine/30)

   			! If an intersection takes place change inter_total from its default
                        ! value of .false. to .true.
                           inter_total = .true.

			! Build wx, wz, unit vector components from point 1 to point 2
                          r1 = polygons(i)%points(j)%r * drfine
                          z1 = polygons(i)%points(j)%z * dzfine
                          r2 = polygons(i)%points(k)%r * drfine
                          z2 = polygons(i)%points(k)%z * dzfine
						
		
                          den = sqrt( (r2-r1)**2 + (z2-z1)**2 )
                          wx = (r2-r1)/den
                          wz = (z2-z1)/den

			! Build gx, gz, unit vector components that point into the metal.
                        ! Because as we move from point 1 to point 2 the metal is on our left
                        ! g = cross(theta,w), where theta is the unit vector that points
                        ! out of the page in the r-z plane.
                           gx = wz
                           gz = -wx

			! Define a point in the cell (xseg,zseg). This seems to work OK
                        ! if the point is on the cell boundary. The logic below relies
                        ! on the fact that the hits... variables are set to -1.d0 if
                        ! no hit occurs on the indicated side of the cell.
      			  if( hitsbottom%r > 0.0 .and. hitsbottom%z > 0.0 ) then
				xseg = hitsbottom%r*drfine
				zseg = hitsbottom%z*dzfine
                           else if( hitstop%r > 0.0 .and. hitstop%z > 0.0 ) then
				xseg = hitstop%r*drfine
				zseg = hitstop%z*dzfine
			   else if( hitsleft%r > 0.0 .and. hitsleft%z > 0.0 ) then
				xseg = hitsleft%r*drfine
				zseg = hitsleft%z*dzfine
			   else if( hitsright%r > 0.0 .and. hitsright%z > 0.0 ) then
				xseg = hitsright%r*drfine
				zseg = hitsright%z*dzfine
			   end if

                        ! set dxstep
                        dxstep=drfine/200.d0

                        ! load each sliced off volume
                        crossings = crossings+1
 			volume(crossings) =  volcell(top*drfine, &
                                 bottom*drfine, left*dzfine, right*dzfine, &
                                 wx, wz, gx, gz, xseg, zseg, dxstep, cellcounter)

!                         if(cellcounter.eq.8101.and.flag.eq.1) then
!
!                tmp=3.14159265359d0*(hitstop%z-left)*dzfine/12.d0*(12.d0*.5* &
!                      (top+bottom)*(top-bottom)*drfine**2+2.d0*(top-bottom)**2*drfine**2)
!!                       pi*dz/12.d0*(12.d0*r0*dr+2.d0*dr**2)
!
!                           write(*,*) 'Volume &&',volume,tmp
!
!                           write(*,*) 'Top: ',hitstop
!                           write(*,*) 'Bottom: ',hitsbottom
!                           write(*,*) 'Left: ',hitsleft
!                           write(*,*) 'Right: ',hitsright
!                           write(*,*) 'drfine, dzfine',drfine,dzfine
!                           write(*,*) 'top,bot,lef,rft',top*drfine, bottom*drfine, left*dzfine, right*dzfine
!                           write(*,*) 'dz outside',hitstop%z-left
!
!
!                         end if


                    end if
                    ! bottom of cell intersection if

                end if
               ! bottom of "last segment if closed polygon" if
                
            end do
!           bottom of polygon segment list

        end do
!   bottom of numpoly list

	! the logical variable "volcorrected" is true if the volume had to
        ! to be corrected, false if not
	volcorrected = inter_total

        ! assign the returned cell volume
        if( inter_total == .false. ) then
            ! the cell is completely inside metal or completely outside metal
            corrected_cell_volume = volume0

        else
            ! at least one shaved-off side occurred; load the accumulated volumes

            ! first, if there were three or more crossings the volume is not going to be
            ! computed correctly
            if(crossings.ge.3) then
              write(*,*) '******************************************************'
              write(*,*) '***********3 or more crossings, bogus volume in cell**',cellcounter
              write(*,*) '******************************************************'
            end if

            ! if the sum of the shaved-off volumes is less than the total volume
            ! then the volume to be included is their sum

            ! if the sum of the shaved-off volumes is greater than the total volume
            ! then the volume to be included is the shaded intersection of the two
            ! shaved-off volumes, which is 2*volume0-volume1-volume2

            volsum=volume(1)+volume(2)

            if(volsum.le.volume0) then
              corrected_cell_volume=volsum
            else
              corrected_cell_volume=2.d0*volume0-volsum
            end if



        end if
    
    return
    end function corrected_cell_volume


!**************************************************************************

    ! returns false if the line segment doesn't intersect with any side of the cell
    ! returns true if it does

     ! variable definitiions:

     ! segpt1:  (point_type)  the first point of the segment (incoming)
     ! segpt2:  (point_type)  the second point of the segment (incoming)

     ! rbottom, rtop, zleft, zright:  The locations (in drfine and
     ! dzfine units) of the edges of the cell  (incoming)

     ! hitsbottom, hitstop, hitsleft, hitsright (point_type): intersection
     ! points in units of drfine and dzfine (outgoing)
     ! If there is no intersection on these particular sides of the volume,
     ! then these have values set to -1.

     ! cellcounter: the cell we're working on, the line number in
     ! collstuff.dat (flag=0) or samplingcells.dat (flag=1)


	logical function doesintersect( segpt1, segpt2, rbottom, &
                                       rtop, zleft, zright, hitsbottom, hitstop, &
                                       hitsleft, hitsright, cellcounter, flag)
	use types_mod
        implicit none

	type(point_type) segpt1, segpt2
	real(8) rbottom, rtop, zleft, zright, s, r, z
        real(8) epsilon,eps,onem,onep ! small floating point limit - adjust to prevent false hits
        real(8) rhits(4),zhits(4),chk

	type(point_type) hitsbottom, hitstop, hitsleft, hitsright
	integer cellcounter,flag,hits

     ! epsilon and the factors used to push lower limits up and upper limits down
     ! if you get zero volume in cells that should have volume, try playing with
     ! this number
        data epsilon/1d-15/

     ! for horizontal and vertical segments use stringent limits to prevent false crossings
     ! for angled segments use slightly generous limits to pick up crossings at cell corners

        if( abs(segpt1%z-segpt2%z)/max(segpt1%z,segpt2%z).lt.epsilon .or. &
            abs(segpt1%r-segpt2%r)/max(segpt1%r,segpt2%r).lt.epsilon ) then
           eps=epsilon
           onem = 1.d0-eps
           onep = 1.d0+eps
        else
           eps=-epsilon
           onem = 1.d0-eps
           onep = 1.d0+eps
        end if


     ! assume at the start that the segment between points 1 and 2 doesn't cross the
     ! cell defined by rbottom, rtop, zleft, and zright.
     ! The logic that follows will correct this if a crossing occurs

       doesintersect = .false.
       hitsbottom%r=-1.d0
       hitsbottom%z=-1.d0
       hitstop%r=-1.d0
       hitstop%z=-1.d0
       hitsright%r=-1.d0
       hitsright%z=-1.d0
       hitsleft%r=-1.d0
       hitsleft%z=-1.d0

     ! set the number of crossings found to zero to start
       hits=0
       rhits(:)=0.d0
       zhits(:)=0.d0

     ! The calculation is based on the parametric form of the line between
     ! the two points segpt1%(r,z) and segpt2%(r,z)
     !  r(s) = (1-s)*r1+s*r2    z(s) = (1-s)*z1 + s*z2  with s between 0 and 1

     ! top of the cell

       if(segpt1%r.ne.segpt2%r) then ! for a top crossing r1=r2 means no crossing

          s=(rtop-segpt1%r)/(segpt2%r-segpt1%r)  ! find the s-value where the line crosses r=rtop
          z=(1.d0-s)*segpt1%z+s*segpt2%z ! find the z-value of the crossing


          ! check for a hit: s between 0 and 1 (on the segment) and z between zleft and zright (in the cell)

          if( s.gt.eps .and. s.lt.onem .and. z.lt.zright*onem .and. z.gt.zleft*onep) then
             hitstop%r=rtop ! crossing at rtop
             hitstop%z=z ! z-value of the crossing at rtop
             doesintersect = .true.
             hits=hits+1
             rhits(hits)=hitstop%r
             zhits(hits)=hitstop%z
!            write(*,*) 'top hit'
          end if

       end if

     ! the other sides of the cell are handled in a similar way

     ! bottom of the cell

       if(segpt1%r.ne.segpt2%r) then

          s=(rbottom-segpt1%r)/(segpt2%r-segpt1%r)
          z=(1.d0-s)*segpt1%z+s*segpt2%z

          ! check for a hit: s between 0 and 1 and z between zleft and zright
          if( s.gt.eps .and. s.lt.onem .and. z.lt.zright*onem .and. z.gt.zleft*onep) then
             hitsbottom%r=rbottom
             hitsbottom%z=z
             doesintersect = .true.
             hits=hits+1
             rhits(hits)=hitsbottom%r
             zhits(hits)=hitsbottom%z
!            write(*,*) 'bottom hit'
          end if

       end if

     ! right side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zright-segpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r


          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.eps .and. s.lt.onem .and. r.lt.rtop*onem .and. r.gt.rbottom*onep) then
             hitsright%r=r
             hitsright%z=zright
             doesintersect = .true.
             hits=hits+1
             rhits(hits)=hitsright%r
             zhits(hits)=hitsright%z
!            write(*,*) 'right hit'
          end if

       end if

     ! left side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zleft-segpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r


          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.eps .and. s.lt.onem .and. r.lt.rtop*onem .and. r.gt.rbottom*onep) then
             hitsleft%r=r
             hitsleft%z=zleft
             doesintersect = .true.
             hits=hits+1
             rhits(hits)=hitsleft%r
             zhits(hits)=hitsleft%z
!            write(*,*) 'left hit'
          end if

       end if

     ! if there are two hits and they are both at the same point, this is
     ! a bogus intersection case, the ones where the cell has an angled
     ! segment ending at one of its corners.

     ! if there is one hit then the segment started on a cell edge -
     ! don't count it as an intersection
         if(hits.eq.1) then

                doesintersect = .false.
                hitsbottom%r=-1.d0
                hitsbottom%z=-1.d0
                hitstop%r=-1.d0
                hitstop%z=-1.d0
                hitsright%r=-1.d0
                hitsright%z=-1.d0
                hitsleft%r=-1.d0
                hitsleft%z=-1.d0

         end if

     ! this is the normal intersection case
         if(hits.eq.2) then
            chk=abs(rhits(1)-rhits(2))+abs(zhits(1)-zhits(1))
            if(chk.lt.epsilon) then

                doesintersect = .false.
                hitsbottom%r=-1.d0
                hitsbottom%z=-1.d0
                hitstop%r=-1.d0
                hitstop%z=-1.d0
                hitsright%r=-1.d0
                hitsright%z=-1.d0
                hitsleft%r=-1.d0
                hitsleft%z=-1.d0
            end if
         end if

     ! If there are three hits, then two of them are identical (at cell corners)
     ! and one needs to be kept and the other eliminated.
     ! If there are four hits, then the segment came in at one
     ! corner and went out at the other. Consolidate these four
     ! hits into just two by eliminating two of them.
     ! Happily, the same logic works for both cases.

         if(hits.eq.3.or.hits.eq.4) then

            ! top=right
            chk=abs(hitstop%r-hitsright%r)+abs(hitstop%z-hitsright%z)
            if(chk.lt.epsilon) then
               hitstop%r=-1.d0
               hitstop%z=-1.d0
            end if

            ! top=left
            chk=abs(hitstop%r-hitsleft%r)+abs(hitstop%z-hitsleft%z)
            if(chk.lt.epsilon) then
               hitstop%r=-1.d0
               hitstop%z=-1.d0
            end if

            ! bottom=left
            chk=abs(hitsbottom%r-hitsleft%r)+abs(hitsbottom%z-hitsleft%z)
            if(chk.lt.epsilon) then
               hitsbottom%r=-1.d0
               hitsbottom%z=-1.d0
            end if
            ! bottom=right

            chk=abs(hitsbottom%r-hitsright%r)+abs(hitsbottom%z-hitsright%z)
            if(chk.lt.epsilon) then
               hitsbottom%r=-1.d0
               hitsbottom%z=-1.d0
            end if

         end if


        return
	end function doesintersect

!**************************************************************************


    ! Volume computation of a cell with part sliced off
    ! find the intersection points 1 and 2 on the edges of the cell
    ! handle horizontal and vertical segments

	! wx, wz -- unit vectors tangent to the segment
	! gx, gz -- unit vectors pointing INTO the metal
	! x, z -- a point on the segment (inside the collision cell)
	! dxstep -- defines tolerance for real equality check (say...drfine/30)

	! cellcounter -- the cell we're on, for testing purposes

        real(8) function volcell(rtop,rbottom,zleft,zright,wx,wz,gx,gz,x,z,dxstep,cellcounter)
        implicit none
        real(8) rtop,rbottom,zleft,zright,wx,wz,gx,gz,x,z,dxstep
	integer cellcounter
        real(8) Volume,drdz,dzdr,chk,r1,r2,z1,z2,dr,dz,r0
        real(8) rmid,rcenter,zcenter,Vtotal,Vright,Vbelow,gdotrc
        real(8) rtmp,z0,zmid,Vcorner

	real(8) :: pi = 3.14159265358979


        ! the crossing segment is horizontal
        if(wx.eq.0.d0) then

          if(gx.gt.0.d0) then
              Volume=pi*(x**2-rbottom**2)*(zright-zleft)
          else
              Volume=pi*(rtop**2-x**2)*(zright-zleft)
          end if

          volcell=Volume
          return
        end if

    	! the crossing segment is vertical
           if(wz.eq.0.d0) then

           if(gz.gt.0.d0) then
              Volume=pi*(rtop**2-rbottom**2)*(z-zleft)
           else
              Volume=pi*(rtop**2-rbottom**2)*(zright-z)
           end if

           volcell=Volume
           return

        end if

	! the crossing segment is neither vertical nor horizontal


	! first, assign to Vtotal the uncorrected volume.
           Vtotal = pi * abs(zright - zleft) * ((rtop)**2 - (rbottom)**2)


        ! find the intersection points 1 and 2 on the edges of the cell
        drdz=wx/wz
        dzdr=wz/wx
        r1=x+drdz*(zleft-z)
        r1=min(r1,rtop)
        r1=max(r1,rbottom)
        r2=x+drdz*(zright-z)
        r2=min(r2,rtop)
        r2=max(r2,rbottom)

        if(r1.gt.r2) then
          rtmp=r1
          r1=r2
          r2=rtmp
        end if

        z1=z+(r1-x)/drdz
        z2=z+(r2-x)/drdz


        ! now find the cell volume by identifying the type of region:
        ! 1. corner or everything but corner
	! 2. quadrilateral caused by vertically cutting down from top to bottom
	! 3. quadrilateral caused by horizontally cutting across from left to right

        ! the test for a corner is: one of r1 or r2 (but not both) lie on
        ! the top or bottom of the cell


        ! defines tolerance for real equality check
        chk=abs(dxstep*1d-7)


        ! if abs(r1-rbottom) < chk & abs(r2-rtop) > chk | abs(r2-rtop) < chk & abs(r1-rbottom)>chk


        ! **** top of corner-vertical-horizontal if block
        if(( abs(r1-rbottom).lt.chk.and.abs(r2-rtop).gt.chk ).or.  &
        (abs(r2-rtop).lt.chk.and.abs(r1-rbottom).gt.chk )) then

            ! it's a corner
            dr=r2-r1
            dz=abs(z2-z1)
            z0=.5d0*dz
            r0=.5d0*(r1+r2)

            ! find a vector that points from the midpoint of the corner
            ! hypotenuse to the cell center
            rmid=.5d0*(r1+r2)
            zmid=.5d0*(z1+z2)
            rcenter=.5d0*(rtop+rbottom)-rmid
            zcenter=.5d0*(zleft+zright)-zmid

            ! find the corner volume
!kluge
!           if(cellcounter.eq.8101) then
!              write(*,*) 'corner: dz,r0,dr',dz,r0,dr
!           end if

            if(rcenter.gt.0.d0) then
                Vcorner=pi*dz/12.d0*(12.d0*r0*dr-2.d0*dr**2)
            else
                Vcorner=pi*dz/12.d0*(12.d0*r0*dr+2.d0*dr**2)
            end if

            gdotrc=rcenter*gx+zcenter*gz

            ! use the dot product between g (points into the metal) and
            ! the vector to the cell center to decide which volume to take

            if(gdotrc.lt.0.d0) then
                Volume=Vtotal-Vcorner
            else
                Volume=Vcorner
            end if

        else if(abs(r1-rbottom).lt.chk.and.abs(r2-rtop).lt.chk) then
            ! the line comes diagonally through the cell from top to bottom


            ! find the right-side volume
            dr=rtop-rbottom
            r0=.5d0*dr+rbottom
            z0=.5d0*(z1+z2)

            Vright=pi*dr/6.d0*(12.d0*r0*(zright-z0)-dzdr*dr**2)

!kluge
!            if(cellcounter.eq.8101) then
!               write(*,*) 'top-bot: dz,r0,dr,dzdr',dz,r0,dr,dzdr
!               write(*,*) 'z1,z2,z0',z1,z2,z0
!            end if

            ! if gz is positive, then we take the volume on the left
            if(gz.gt.0.d0) then
                Volume=Vtotal-Vright
            else
                Volume=Vright
            end if

        else
            ! line comes diagonally across the cell from left to right
	
            ! find the lower volume
            dz=zright-zleft
            z0=.5d0*dz
            r0=.5d0*(r1+r2)

            Vbelow=pi/12.d0*dz*(drdz**2*dz**2+12.d0*(r0**2-rbottom**2))

!kluge
!            if(cellcounter.eq.8101) then
!               write(*,*) 'lef-rgt: dz,r0,drdz',dz,r0,dr,dzdr
!            end if
!
		! if gx is negative, then we take the volume above the line
		if(gx.lt.0.d0) then
                Volume=Vtotal-Vbelow
            else
                Volume=Vbelow
            end if


        end if
        ! **** bottom of corner-vertical-horizontal if block

        volcell=Volume
      
    return
    end function volcell

end module volume_calc_mod


!**** fgvolumecalc.f90 ****
