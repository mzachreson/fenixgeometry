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
        real(8) bottom, top, left, right
	logical volcorrected
	integer cellcounter, flag

        real(8) volume

	! vector variables:
	real(8) wx, wz, gx, gz
		
	! used to pass in the point on the segment into "volcell"
	real(8) onsegr, onsegz

	! temporary variables for calculation of vectors
	real(8) r1, r2, z1, z2
	real(8) den

	integer i, j, k
	logical inter_thissegment, inter_thissegment2
	logical inter_total

	real(8) shouldbe ! temperary for testing

        type(point_type) hitsbottom, hitstop, hitsleft, hitsright
		


        ! I want to create "integral_regions", which specify the regions which I want to do in integral
        ! just like in Math 113, to find the volume of the area revolved around the axis...


        ! step 1: loop over all the line segments and find all the places where there's an intersection
		inter_total = .false.
        ! start the corrected volume out at 0. Everytime a segment cuts through the cell and
        ! we calculate the volume on the outside of the cut, add it to volume. This allows
        ! for a proper calculation when two segments cut through a cell

        do i = 1,npoly
            do j = 1,polygons(i)%numpoints
                if( j == polygons(i)%numpoints ) then
                    k = 1
                else
                    k = j + 1
                end if


                ! only do the last one if it's a closed polygon
                if( j < polygons(i)%numpoints .or. polygons(i)%polytype == 1 ) then

                    inter_thissegment2= doesintersect( polygons(i)%points(j), &
                                                polygons(i)%points(k), &
                                                bottom, top, left, right, &
                                                hitsbottom, hitstop, hitsleft, hitsright, cellcounter, flag)
! test my version
                    inter_thissegment = doesintersect2( polygons(i)%points(j), &
                                                polygons(i)%points(k), &
                                                bottom, top, left, right, &
                                                hitsbottom, hitstop, hitsleft, hitsright, cellcounter, flag)

                     if(inter_thissegment.ne.inter_thissegment2) then
                       write(*,*) '$$ ',cellcounter,inter_thissegment,inter_thissegment2,i,j,k
                     end if

                     if(cellcounter.eq.21501.and.flag.eq.1) then

                       write(*,*) '## ',cellcounter,inter_thissegment,inter_thissegment2,i,j,k
                     end if


! test my  version

                    if( inter_thissegment == .true. ) then
                        ! then calculate a special corrected volume with "volcell"

			! wx, wz -- unit vectors tangent to the segment
			! gx, gz -- unit vectors pointing INTO the metal
			! x, z -- a point on the segment (inside the collision cell)
			! dxstep -- defines tolerance for real equality check (say...drfine/30)

   			! this is so we know that AT LEAST ONE of the segments intersected.
				inter_total = .true.

			! defining wx, wz
				r1 = polygons(i)%points(j)%r * drfine
				z1 = polygons(i)%points(j)%z * dzfine
				r2 = polygons(i)%points(k)%r * drfine
				z2 = polygons(i)%points(k)%z * dzfine
						
		
				den = sqrt( (r2-r1)**2 + (z2-z1)**2 )
				wx = (r2-r1)/den
				wz = (z2-z1)/den

			!defining gx, gz, so that it's perpendicular to <wx, wz> and that it points
			! toward the INTERIOR of the metal, which is on the left as you follow the points around.
				gx = wz
				gz = -wx

			! finding onsegr, onsegz
			! 	here I am assuming it's okay that the included point is on the boundary
			!   (and not in between, inside the cell.) I hope that's right.
      			  if( hitsbottom%r > 0.0 .and. hitsbottom%z > 0.0 ) then
				onsegr = hitsbottom%r*drfine
				onsegz = hitsbottom%z*dzfine
                           else if( hitstop%r > 0.0 .and. hitstop%z > 0.0 ) then
				onsegr = hitstop%r*drfine
				onsegz = hitstop%z*dzfine
			   else if( hitsleft%r > 0.0 .and. hitsleft%z > 0.0 ) then
				onsegr = hitsleft%r*drfine
				onsegz = hitsleft%z*dzfine
			   else if( hitsright%r > 0.0 .and. hitsright%z > 0.0 ) then
				onsegr = hitsright%r*drfine
				onsegz = hitsright%z*dzfine
			   end if

                        ! add up all of the sliced off volumes
 			volume = volume + volcell(top*drfine, bottom*drfine, left*dzfine, right*dzfine, &
                                 wx, wz, gx, gz, onsegr, onsegz, drfine/50.0, cellcounter)
                        if(cellcounter.eq.8101.and.flag.eq.1) write(*,*) 'Volume&&',volume


                    end if
                end if

                
            end do

        end do

		! for testing purposes, this passed in variable is set
		volcorrected = inter_total

        ! assign the return value
        if( inter_total == .false. ) then
            ! then either the box is completely inside a metal or it's completely outside a metal
			! for now, we just assume completely outside a metal
            corrected_cell_volume = pi * abs(right*dzfine - left*dzfine) * ((top*drfine)**2 - (bottom*drfine)**2)

        else
            corrected_cell_volume = volume

        end if
    
    return
    end function corrected_cell_volume





    ! returns false if the line segment doesn't intersect with any side of the cell
    ! returns true if it does

    ! variable definitiions:

    ! segpt1:  (point_type)  the first point of the segment (incoming)
    ! segpt2:  (point_type)  the second point of the segment (incoming)
    ! bottom, top, left, right:  The locations (in dr/dzfine units) of the edges of the cell  (incoming)
    ! hitsbottom, hitstop, hitsleft, hitsrigth:  (OUTGOING)
    !      the point_type points where the intersection occurs.
    !      If there is no intersection on these particular sides of the volume, then these have negative values
	! cellcounter: the cell we're on (the line number in "collcells.dat")

	! NOTE: this whole function uses dz/drfine units

	logical function doesintersect( segpt1, segpt2, bottom, top, &
                   left, right, hitsbottom, hitstop, hitsleft, hitsright, cellcounter, flag)
		use types_mod
	    implicit none

		type(point_type) segpt1, segpt2
		real(8) bottom, top, left, right
		real(8) rlesser, rgreater, zlesser, zgreater
		type(point_type) hitsbottom, hitstop, hitsleft, hitsright
		integer cellcounter, flag
        logical bool_bottom, bool_top, bool_left, bool_right
		
		real(8) m  ! slope of the line
		
    	! equations:
    	! box1:  r = bottom
    	! box2:  r = top
    	! box3:  z = left
    	! box4:  z = right
    	! segment:  m = ((segpt2%r-segpt1%r)/(segpt2%z-segpt1%z))
    	!           r = m(z - z0) + r0

    	!special case 1:  slope is 0
    	if( abs(segpt2%r-segpt1%r) < 1d-20) then


		    !finding out which points are lesser, and which ones are greater
            if( segpt1%z < segpt2%z ) then
			    zlesser = segpt1%z
			    zgreater = segpt2%z
		    else
			    zlesser = segpt2%z
			    zgreater = segpt1%z
		    end if

            ! with a horizontle line, it's impossible to intersect with the bottom or top
            ! in a way that will affect the volume.
            hitsbottom%r = -1.0
            hitsbottom%z = -1.0
            hitstop%r = -1.0
            hitstop%r = -1.0
            bool_bottom = .false.
            bool_top = .false.

            ! but it can hit the left and right:
            hitsleft%r = dble(segpt1%r)
            hitsleft%z = dble(left)

    		hitsright%r = dble(segpt1%r)
            hitsright%z = dble(right)

			! set these to true initially, and then check if it's false
            bool_left = .true.
            bool_right = .true.

    		if(  .not. ( hitsleft%r < top .and. hitsleft%r > bottom  &   
                 .and.   hitsleft%z < zgreater .and. hitsleft%z >= zlesser ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsleft%r = -1.0
                hitsleft%z = -1.0
                bool_left = .false.
                        
            end if

    		if(  .not. ( hitsright%r < top      .and. hitsright%r > bottom  &   
                 .and.   hitsright%z <= zgreater .and. hitsright%z > zlesser ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsright%r = -1.0
                hitsright%z = -1.0
                bool_right = .false.
                        
            end if
    	
    	








    	!special case 2:  slope is infinity, line is vertical
    	else if( abs(segpt2%z-segpt1%z) < 1d-10 ) then


	    !finding out which points are lesser, and which ones are greater
	    if( segpt1%r < segpt2%r ) then
		    rlesser = segpt1%r
		    rgreater = segpt2%r
	    else
		    rlesser = segpt2%r
		    rgreater = segpt1%r
	    end if

            !with a virtical line, it's impossible for the left and right to be intersected
            ! in a way that affects the volume
            hitsleft%r = -1.0
            hitsleft%z = -1.0
            hitsright%r = -1.0
            hitsright%r = -1.0
            bool_left = .false.
            bool_right = .false.

            ! but it can hit the top and bottom
            hitstop%r = dble(top)
            hitstop%z = dble(segpt1%z)

     	    hitsbottom%r = dble(bottom)
            hitsbottom%z = dble(segpt1%z)
    	
			! set these to true initially, and then check if it's false
            bool_bottom = .true.
            bool_top = .true.

    		if(  .not. ( hitsbottom%r < rgreater .and. hitsbottom%r >= rlesser &  
    			 .and.   hitsbottom%z < right    .and. hitsbottom%z > left    ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsbottom%r = -1.0
                hitsbottom%z = -1.0
                bool_bottom = .false.
                        
            end if

    		if(  .not. ( hitstop%r <= rgreater .and. hitstop%r > rlesser &   
    			 .and.   hitstop%z < right    .and. hitstop%z > left    ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitstop%r = -1.0
                hitstop%z = -1.0
                bool_top = .false.
                        
            end if








    	
    	! slope is a finite number 0 < m < Infinity
    	else


		    !finding out which points are lesser, and which ones are greater
		    if( segpt1%r < segpt2%r ) then
			    rlesser = segpt1%r
			    rgreater = segpt2%r
		    else
			    rlesser = segpt2%r
			    rgreater = segpt1%r
		    end if
		
		
		    if( segpt1%z < segpt2%z ) then
			    zlesser = segpt1%z
			    zgreater = segpt2%z
		    else
			    zlesser = segpt2%z
			    zgreater = segpt1%z
		    end if


    		m = (dble(segpt2%r-segpt1%r)/dble(segpt2%z-segpt1%z));
    	
    		! finding the z or r value of intersection
    		hitsbottom%z = dble(bottom - segpt1%r) / m + dble(segpt1%z)
            hitsbottom%r = m*(dble(hitsbottom%z - segpt1%z)) + dble(segpt1%r)

    		hitstop%z = dble(top - segpt1%r) / m + dble(segpt1%z)
            hitstop%r = m*(dble(hitstop%z - segpt1%z)) + dble(segpt1%r)
    	
    		hitsleft%r = m*(dble(left-segpt1%z)) + dble(segpt1%r)
            hitsleft%z = dble(hitsleft%r - segpt1%r) / m + dble(segpt1%z)

    		hitsright%r = m*(dble(right-segpt1%z)) + dble(segpt1%r)
            hitsright%z = dble(hitsright%r - segpt1%r) / m + dble(segpt1%z)


			! set these to true initially, and then check if it's false
            bool_bottom = .true.
            bool_top = .true.
            bool_left = .true.
            bool_right = .true.

    		if(  .not. ( hitsbottom%r < (rgreater+1d-10) .and. hitsbottom%r > (rlesser-1d-10) &  
                 .and.   hitsbottom%z < (zgreater+1d-10) .and. hitsbottom%z > (zlesser-1d-10) &
    			 .and.   hitsbottom%z < (right+1d-10)    .and. hitsbottom%z > (left-1d-10)    ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsbottom%r = -1.0
                hitsbottom%z = -1.0
                bool_bottom = .false.
                        
            end if

    		if(  .not. ( hitstop%r < (rgreater+1d-10) .and. hitstop%r > (rlesser-1d-10) &   
                 .and.   hitstop%z < (zgreater+1d-10) .and. hitstop%z > (zlesser-1d-10) &
    			 .and.   hitstop%z < (right+1d-10)    .and. hitstop%z > (left-1d-10)    ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitstop%r = -1.0
                hitstop%z = -1.0
                bool_top = .false.
                        
            end if

    		if(  .not. ( hitsleft%r < (rgreater+1d-10) .and. hitsleft%r > (rlesser-1d-10) &
                 .and.   hitsleft%r < (top+1d-10)      .and. hitsleft%r > (bottom-1d-10)  &   
                 .and.   hitsleft%z < (zgreater+1d-10) .and. hitsleft%z > (zlesser-1d-10) ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsleft%r = -1.0
                hitsleft%z = -1.0
                bool_left = .false.
                        
            end if

    		if(  .not. ( hitsright%r < (rgreater+1d-10) .and. hitsright%r > (rlesser-1d-10) &
                 .and.   hitsright%r < (top+1d-10)      .and. hitsright%r > (bottom-1d-10)  &   
                 .and.   hitsright%z < (zgreater+1d-10) .and. hitsright%z > (zlesser-1d-10) ) ) then

                !if it's not in the span of the segment or the span of the box, then it doesn't hit the box
                hitsright%r = -1.0
                hitsright%z = -1.0
                bool_right = .false.
                        
            end if
 	
    	
    	end if


        ! setting the return value
        if(       bool_bottom == .false. .and. bool_top == .false. &
            .and. bool_left == .false. .and. bool_right == .false. ) then

            doesintersect = .false.
        else
            doesintersect = .true.
        end if


	end function doesintersect




    ! returns false if the line segment doesn't intersect with any side of the cell
    ! returns true if it does

	! variable definitiions:

	! segpt1:  (point_type)  the first point of the segment (incoming)
	! segpt2:  (point_type)  the second point of the segment (incoming)
	! rbottom, rtop, zleft, zright:  The locations (in drfine and dzfine units) of the edges of the cell  (incoming)
    ! hitsbottom, hitstop, hitsleft, hitsright:  (OUTGOING) in units of drfine and dzfine
    !      the point_type points where the intersection occurs.
    !      If there is no intersection on these particular sides of the volume, then these have values of -1.
	! cellcounter: the cell we're on (the line number in "collcells.dat")


	logical function doesintersect2( segpt1, segpt2, rbottom, &
           rtop, zleft, zright, hitsbottom, hitstop, hitsleft, hitsright, cellcounter, flag)
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

        if(cellcounter.eq.21501.and.flag.eq.1) then
            write(*,*) '%%% r1,r2',segpt1%r,segpt2%r,segpt1%r-segpt2%r
            write(*,*) 'z1,z2: ',segpt1%z,segpt2%z,segpt1%z-segpt2%z
            write(*,*) 'eps',eps
        end if

		

       ! assume at the start that the segment between points 1 and 2 doesn't cross the cell defined by rbottom, rtop, zleft, and zright.
       ! The logic below will correct this if a crossing occurs

       doesintersect2 = .false.
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

       ! The calculation is based on the parametric form of the line between the two points segpt1%(r,z) and segpt2%(r,z)
       !  r(s) = (1-s)*r1+s*r2    z(s) = (1-s)*z1 + s*z2  with s between 0 and 1

       ! top of the cell

       if(segpt1%r.ne.segpt2%r) then ! for a top crossing r1=r2 means no crossing

          s=(rtop-segpt1%r)/(segpt2%r-segpt1%r)  ! find the s-value where the line crosses r=rtop
          z=(1.d0-s)*segpt1%z+s*segpt2%z ! find the z-value of the crossing


          ! check for a hit: s between 0 and 1 (on the segment) and z between zleft and zright (in the cell)

          if( s.gt.eps .and. s.lt.onem .and. z.lt.zright*onem .and. z.gt.zleft*onep) then
             hitstop%r=rtop ! crossing at rtop
             hitstop%z=z ! z-value of the crossing at rtop
             doesintersect2 = .true.
             hits=hits+1
             rhits(hits)=hitstop%r
             zhits(hits)=hitstop%z
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
             doesintersect2 = .true.
             hits=hits+1
             rhits(hits)=hitsbottom%r
             zhits(hits)=hitsbottom%z
          end if

       end if

       ! right side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zright-segpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r


          if(cellcounter.eq.21501.and.flag.eq.1) then
            write(*,*) 's,r,zright',s,r,zright
            write(*,*) 'eps',eps

          end if

          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.eps .and. s.lt.onem .and. r.lt.rtop*onem .and. r.gt.rbottom*onep) then
             hitsright%r=r
             hitsright%z=zright
             doesintersect2 = .true.
             hits=hits+1
             rhits(hits)=hitsright%r
             zhits(hits)=hitsright%z
          end if

       end if

       ! left side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zleft-segpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r

          if(cellcounter.eq.21501.and.flag.eq.1) then
!           write(*,*) 's-onem,r-rbottom*onep',s-onem,r-rbottom*onep

!           write(*,*) 'eps',eps
          end if

          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.eps .and. s.lt.onem .and. r.lt.rtop*onem .and. r.gt.rbottom*onep) then
             hitsleft%r=r
             hitsleft%z=zleft
             doesintersect2 = .true.
             hits=hits+1
             rhits(hits)=hitsleft%r
             zhits(hits)=hitsleft%z
          end if

       end if

! if there are two hits and they are both at the same point, this is
! a bogus intersection case, the ones where the cell has an angled
! segment ending at one of its corners.

         if(hits.eq.2) then
            chk=abs(rhits(1)-rhits(2))+abs(zhits(1)-zhits(1))
            if(chk.lt.epsilon) then

                doesintersect2 = .false.
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

! if there are three hits, then two of them are identical (at cell corners)
! and one needs to be kept and the other eliminated

         if(hits.eq.3) then

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


              if(cellcounter.eq.21501.and.flag.eq.1) then
               write(*,*) 'Check, 21501 '
               write(*,*) 'r1,z1: ',segpt1%r,segpt1%z
               write(*,*) 'r2,z2: ',segpt2%r,segpt2%z
               write(*,*) 'hit top (r,z)',hitstop%r,hitstop%z
               write(*,*) 'hit bottom (r,z)',hitsbottom%r,hitsbottom%z
               write(*,*) 'hit right (r,z)',hitsright%r,hitsright%z
               write(*,*) 'hit left (r,z)',hitsleft%r,hitsleft%z
              write(*,*) 'doesintersect2', doesintersect2
              end if

        return
	end function doesintersect2





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
            z0=.5*dz
            r0=.5*(r1+r2)

            ! find a vector that points from the midpoint of the corner
            ! hypotenuse to the cell center
            rmid=.5*(r1+r2)
            zmid=.5*(z1+z2)
            rcenter=.5*(rtop+rbottom)-rmid
            zcenter=.5*(zleft+zright)-zmid



            ! find the corner volume

            if(rcenter.gt.0.d0) then
                Vcorner=pi*dz/12*(12*r0*dr-2*dr**2)
            else
                Vcorner=pi*dz/12*(12*r0*dr+2*dr**2)
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
            r0=.5*dr+rbottom
            z0=.5*(z1+z2)

            Vright=pi*dr/6.d0*(12.d0*r0*(zright-z0)-dzdr*dr**2)

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
            z0=.5*dz
            r0=.5*(r1+r2)

            Vbelow=pi/12.d0*dz*(drdz**2*dz**2+12.d0*(r0**2-rbottom**2))

			! if gx is negative, then we take the volume above the line
			if(gx.lt.0.d0) then
                Volume=Vtotal-Vbelow
            else
                Volume=Vbelow
            end if


        end if
        ! **** bottom of corner-vertical-horizontal if block

        volcell=Volume
      
    end function volcell







    !---cell_volume -- may be modified in the future to
    !               -- account for part of the cell inside
    !               -- a metal

    !               -- Assumes a complete circle (2 pi radians)
    !               -- Assumes values in units of drfine/dzfine
    real(8) function cell_volume( bottom, top, left, right, drfine, dzfine)

        implicit none
        real(8) bottom, top, left, right, drfine, dzfine
        real(8) :: pi = 3.14159265358979

        cell_volume = pi * abs(right*dzfine - left*dzfine) * ((top*drfine)**2 - (bottom*drfine)**2)

    end function cell_volume


! Currently, this function isn't being used:

    ! given a point, and a box, returns whether the point is inside the box (not including the edge)
!    logical function pointinside( pr, pz, bottom, top, left, right )
!        implicit none

!        real(8) pr, pz, bottom, top, left, right

!        if( pr < top .and. pr > bottom .and. pz < right .and. pz > left ) then
!            pointinside = .true.
!        else
!            pointinside = .false.
!        end if

!    end function pointinside


end module volume_calc_mod


