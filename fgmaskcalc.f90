!*** fgmaskcalc.f90 ***





module mask_calc_mod
	implicit none

	contains



!This function is called for each polygon


! The following are now encapsulated in the polygon_type variable "polygon":
! 		INCOMING: zpoly:  Array of z-points in the polygon (in meters)
! 		INCOMING: rpoly:  Array of r-points in the polygon (in meters)
! 		INCOMING: polytype: polytype of this polygon.
! 		INCOMING: npoints:  number of points in the polygon
! 		INCOMING: reflection: the reflection code for the polygon
! 		INCOMING: temp:  the temperature of the polygon

! INCOMING: step:  the width of the painted region in meters. step = 8*vth*tau
! INCOMING: iwall: initialized to zero before any polygons are processed, but as the code
		! moves from polygon to polygon it must be iterated in steps of 2 to make EVEN (for a wall) numbered mask codes
! INCOMING: ijoint: initialized to -1, before any polygons are processed, but as the code
		! moves from polygon to polygon it must be iterated in steps of 2 to make ODD (for a joint) numbered mask codes 
! INCOMING: imask: initialized to zero before any polygons are processed, but as the code
		! moves from polygon to polygon it must be iterated in steps of 1, to count mask codes.

! The following are now encapsulated in the mask_type variable "masks" (located in "data_mod")
! 		OUTGOING: w1, g1, w2, g2:  2D arrays that specify tangent (w1 and w2) and perpendicular (g1 and g2)
					! unit vectors for each mask code.
! 		OUTGOING: s0:   2D array, specifies the end point position of the mask code (r and z) (in meters)
! 		OUTGOING: reflectmask:  the reflection code transfered to the mask code
! 		OUTGOING: tempmask the temperature transferred to the mask code.

! The "findcell" array is located in "data_mod" and is already created by the time this subroutine is called.


	subroutine maskpaint( polygon, step, iwall, ijoint, imask )

		use data_mod
		use types_mod
        implicit none

		! passed in values
		type(polygon_type) polygon
        real(8) step
		integer iwall, ijoint, imask

		! used to be passed in (but now info is contained in "polygon"):
		integer npoints
        real(8) zpoly(polygon%numpoints),rpoly(polygon%numpoints)
        integer polytype


		! NOTE: for these "w1, g1", etc, values, 1 is Z, 2 is R
        real(8) w1(totnumpoints,2),g1(totnumpoints,2),w2(totnumpoints,2),g2(totnumpoints,2)
		real(8) s0(totnumpoints,2)
        integer mcode(totnumpoints) ! the mask codes

		real(8) reflection, temp
        real(8) reflectmask(totnumpoints),tempmask(totnumpoints)

		! local values
        real(8) theta,rho,theta0,drho,dtheta,cs,sn
        real(8) rpt,zpt,d,wx,wz,gx,gz,s0z,s0x
        integer i,j,k,m,n,ii,jj,kk
        integer numsteps,ntheta,maskcode,nline
		real(8) :: pi = 3.14159265358979

		integer imaskstart, imaskend

		imaskstart = imask + 1



		!adapting the inputted values to the OLD inputted values
		npoints = polygon%numpoints
		do i = 1,npoints
			rpoly(i) = polygon%points(i)%r * drfine
			zpoly(i) = polygon%points(i)%z * dzfine
		end do
		polytype = polygon%polytype
		temp = polygon%temperature
		reflection = polygon%reflect

		! initializing these arrays so they're zero if not used:
		do i = 1, totnumpoints
			w1(i,1) = 0.0
			g1(i,1) = 0.0
			w2(i,1) = 0.0
			g2(i,1) = 0.0
			s0(i,1) = 0.0

			w1(i,2) = 0.0
			g1(i,2) = 0.0
			w2(i,2) = 0.0
			g2(i,2) = 0.0
			s0(i,2) = 0.0

			reflectmask(i) = 0.0
			tempmask(i) = 0.0
		end do

		! NOW WE START
        numsteps=40
        drho=step/real(numsteps)
        ! set ntheta to be about 2*pi*numsteps
        ntheta=numsteps*6
        dtheta=2.*pi/real(ntheta)
        theta0=.5*dtheta


        ! run over the polygon and assign some mask codes


        ! process each point as follows:

        ! If polytype=0 then the polygon is open, meaning that it connects
        ! to the computation boundary. In this case we make a wall mask from
        ! the boundary connection to halfway down the line to the next vertex,
        ! at the beginning and at the end. We also put a circle of mask codes
        ! centered at each vertex. The polygon must start and end at
        ! the computation boundary and the order of points must be such that
        ! as you travel along it from the first point to the last the interior
        ! of the computation (where particles are) is on your right.

        ! If polytype=1 then the polygon is closed and each vertex is processed
        ! as a joint with the joint extending from the vertex to the halfway
        ! point to each neighboring vertex. A circle of mask codes is also
        ! painted centered at each vertex. Note that for this polytype the
        ! last point is assumed to be the same as the first one, so don't
        ! enter the last point in the input list.

        ! This logic is implemented by having an if-block to separately
        ! handle the first and last points for both polytypes. The rest
        ! of the points are then handled the same in both cases.


        if(polytype.eq.0) then
        

            ! open polygon case - handle the beginning and the end specially here, then
            ! handle the rest of the points in the common polytype section

            ! process the first point with a wall section extending
            ! toward the next point
            n=1
            ! increment the mask code counter
            imask=imask+1
            ! increment the mask code itself, even numbers
            ! mean walls
            iwall=iwall+2
            maskcode=iwall
            mcode(imask)=maskcode


            ! load w1 and g1 for this segment then paint the segment

            w1(imask,1)=zpoly(2)-zpoly(1)
            w1(imask,2)=rpoly(2)-rpoly(1)
            ! normalize it
            d=sqrt(w1(imask,1)**2+w1(imask,2)**2)
            w1(imask,1)=w1(imask,1)/d
            w1(imask,2)=w1(imask,2)/d

            !  g points into the metal, g=cross(yhat,w)
            g1(imask,1)=-w1(imask,2)
            g1(imask,2)=w1(imask,1)

            ! w2 and g2 must be zero for a wall segment
            w2(imask,1)=0.
            w2(imask,2)=0.
            g2(imask,1)=0.
            g2(imask,2)=0.

            ! load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w1(imask,1)
            wx=w1(imask,2)
            gz=g1(imask,1)
            gx=g1(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            ! load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                wx,wz,gx,gz,maskcode,s0z,s0x)

            !c process the last point with a wall section extending
            !c toward the point behind it
            n=npoints
            !c increment the mask code counter
            imask=imask+1
            !c increment the mask code itself, even numbers
            !c mean walls
            iwall=iwall+2
            maskcode=iwall
            mcode(imask)=maskcode


            !c load w1 and g1 for this segment then paint the segment

            w1(imask,1)=zpoly(n-1)-zpoly(n)
            w1(imask,2)=rpoly(n-1)-rpoly(n)
            !c normalize it
            d=sqrt(w1(imask,1)**2+w1(imask,2)**2)
            w1(imask,1)=w1(imask,1)/d
            w1(imask,2)=w1(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w),
            !c  but I am doing this segment backwards from the read-in order
            !c  so the signs are flipped
            g1(imask,1)= w1(imask,2)
            g1(imask,2)=-w1(imask,1)

            !c w2 and g2 must be zero for a wall segment
            w2(imask,1)=0.
            w2(imask,2)=0.
            g2(imask,1)=0.
            g2(imask,2)=0.


            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w1(imask,1)
            wx=w1(imask,2)
            gz=g1(imask,1)
            gx=g1(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                        wx,wz,gx,gz,maskcode, & 
                        s0z,s0x)



        else if(polytype.eq.1) then
            !c closed polygon case


            !c process the first point with a joint section extending
            !c toward the next point (w2 and g2)  and toward the
            !c last point (w1 and g1)
            n=1
            !c increment the mask code counter
            imask=imask+1
            !c increment the mask code itself, odd numbers
            !c mean joints
            ijoint=ijoint+2
            maskcode=ijoint
            mcode(imask)=maskcode

            !c ********* 1->2 ************ (w2 and g2)
            !c load w2 and g2 for this segment then paint the segment

            w2(imask,1)=zpoly(2)-zpoly(1)
            w2(imask,2)=rpoly(2)-rpoly(1)
            !c normalize it
            d=sqrt(w2(imask,1)**2+w2(imask,2)**2)
            w2(imask,1)=w2(imask,1)/d
            w2(imask,2)=w2(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w)
            g2(imask,1)=-w2(imask,2)
            g2(imask,2)=w2(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w2(imask,1)
            wx=w2(imask,2)
            gz=g2(imask,1)
            gx=g2(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0, &
                        wx,wz,gx,gz,maskcode, &
                        s0z,s0x)


            !c ********* 1->npoints ************ (w1 and g1)
            !c load w1 and g1 for this segment then paint the segment

            w1(imask,1)=zpoly(npoints)-zpoly(1)
            w1(imask,2)=rpoly(npoints)-rpoly(1)
            !c normalize it
            d=sqrt(w1(imask,1)**2+w1(imask,2)**2)
            w1(imask,1)=w1(imask,1)/d
            w1(imask,2)=w1(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w), but
            !c we are running backwards from 1->npoints, so flip the sign
            g1(imask,1)= w1(imask,2)
            g1(imask,2)=-w1(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w1(imask,1)
            wx=w1(imask,2)
            gz=g1(imask,1)
            gx=g1(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0, &
                wx,wz,gx,gz,maskcode, &
                s0z,s0x)



            !c process the last point with a joint section extending
            !c toward the first point and toward the next-to-last point
            n=npoints
            !c increment the mask code counter
            imask=imask+1
            !c increment the mask code itself, odd numbers
            !c mean joints
            ijoint=ijoint+2
            maskcode=ijoint
            mcode(imask)=maskcode

            !c ********* npoints->1 ****** (w2 and g2)
            !c load w2 and g2 for this segment then paint the segment

            w2(imask,1)=zpoly(1)-zpoly(npoints)
            w2(imask,2)=rpoly(1)-rpoly(npoints)
            !c normalize it
            d=sqrt(w2(imask,1)**2+w2(imask,2)**2)
            w2(imask,1)=w2(imask,1)/d
            w2(imask,2)=w2(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w)
            g2(imask,1)=-w2(imask,2)
            g2(imask,2)=w2(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w2(imask,1)
            wx=w2(imask,2)
            gz=g2(imask,1)
            gx=g2(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                wx,wz,gx,gz,maskcode, &
                s0z,s0x)


            !c ********* npoints-> npoints-1 *** (w1 and g1)
            !c load w1 and g1 for this segment then paint the segment

            w1(imask,1)=zpoly(npoints-1)-zpoly(npoints)
            w1(imask,2)=rpoly(npoints-1)-rpoly(npoints)
            !c normalize it
            d=sqrt(w1(imask,1)**2+w1(imask,2)**2)
            w1(imask,1)=w1(imask,1)/d
            w1(imask,2)=w1(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w), but
            !c we are running backwards from npoints->npoints-1, so flip the sign
            g1(imask,1)= w1(imask,2)
            g1(imask,2)=-w1(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w1(imask,1)
            wx=w1(imask,2)
            gz=g1(imask,1)
            gx=g1(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                wx,wz,gx,gz,maskcode, &
                s0z,s0x)



        else
            !c illegal polytype
            write(*,*) 'Illegal polytype - ',polytype
            stop 444

        end if


        !c now that the first and last points have been specially handled
        !c we can do the other points with a loop, all done as joints

        !c top of remaining points loop
        do n=2,npoints-1


            !c increment the mask code counter
            imask=imask+1
            !c increment the mask code itself, odd numbers
            !c mean joints
            ijoint=ijoint+2
            maskcode=ijoint
            mcode(imask)=maskcode

            !c ********* n->n+1 ******
            !c load w2 and g2 for this segment then paint the segment

            w2(imask,1)=zpoly(n+1)-zpoly(n)
            w2(imask,2)=rpoly(n+1)-rpoly(n)
            !c normalize it
            d=sqrt(w2(imask,1)**2+w2(imask,2)**2)
            w2(imask,1)=w2(imask,1)/d
            w2(imask,2)=w2(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w)
            g2(imask,1)=-w2(imask,2)
            g2(imask,2)=w2(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w2(imask,1)
            wx=w2(imask,2)
            gz=g2(imask,1)
            gx=g2(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                        wx,wz,gx,gz,maskcode,  &
                        s0z,s0x)


            !c ********* n->n-1 ******
            !c load w1 and g1 for this segment then paint the segment

            w1(imask,1)=zpoly(n-1)-zpoly(n)
            w1(imask,2)=rpoly(n-1)-rpoly(n)
            !c normalize it
            d=sqrt(w1(imask,1)**2+w1(imask,2)**2)
            w1(imask,1)=w1(imask,1)/d
            w1(imask,2)=w1(imask,2)/d

            !c  g points into the metal, g=cross(yhat,w), but
            !c we are running backwards from npoints->npoints-1, so flip the sign
            g1(imask,1)= w1(imask,2)
            g1(imask,2)=-w1(imask,1)

            !c load s0
            s0(imask,1)=zpoly(n)
            s0(imask,2)=rpoly(n)

            wz=w1(imask,1)
            wx=w1(imask,2)
            gz=g1(imask,1)
            gx=g1(imask,2)
            s0z=s0(imask,1)
            s0x=s0(imask,2)

            !c load the reflection coefficient and the temperature for this mask code
            reflectmask(imask)=reflection
            tempmask(imask)=temp


            call painter(d,drho,numsteps,ntheta,dtheta,theta0,  &
                        wx,wz,gx,gz,maskcode, &
                        s0z,s0x)

        end do
        ! end of remaining points loop




		imaskend = imask

		! before we finish, we convert the formerly-returned values of 
		! "wx, wz, gx, gz, s0 masktemp, maskreflect, mcode"
		! to what is now returned, AKA, the "mask_type" array "masks"

		do i = imaskstart,imaskend

			masks(i)%code = mcode(i)
			masks(i)%reflect = reflectmask(i)
			masks(i)%s0r = s0(i,2)
			masks(i)%s0z = s0(i,1)
			masks(i)%w1r = w1(i,2)
			masks(i)%w1z = w1(i,1)
			masks(i)%g1r = g1(i,2)
			masks(i)%g1z = g1(i,1)
			masks(i)%w2r = w2(i,2)
			masks(i)%w2z = w2(i,1)
			masks(i)%g2r = g2(i,2)
			masks(i)%g2z = g2(i,1)
			masks(i)%temperature = tempmask(i)

			if( mcode(i) > maxmaskcode ) then
				maxmaskcode = mcode(i)
			end if

		end do
		

    end subroutine maskpaint


!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************
!*********************************


    subroutine painter(d,drho,numsteps,ntheta,dtheta,theta0, &
                        wx,wz,gx,gz,maskcode,z0,r0)
		use data_mod
        implicit none

        real(8) d,drho,wx,wz,gx,gz,dtheta,cs,sn,rpt,zpt,rho,ds
        real(8) theta,theta0,z0,r0

        integer ii,jj,numsteps,ntheta,nline,i,j
        integer maskcode
		integer q,r,s



        !c draw the circle
        do ii=1,ntheta

            theta=theta0+(ii-1)*dtheta
            cs=cos(theta)
            sn=sin(theta)
		
            do jj=1,numsteps
                rho=.5*drho+(jj-1)*drho
                zpt=z0+rho*cs
                rpt=r0+rho*sn

                ! top of in-bounds if
                if(rpt.lt.(rmax*drfine).and.rpt.gt.(rmin*drfine).and. &
                    zpt.gt.(zmin*dzfine).and.zpt.lt.(zmax*dzfine)) then

                    i = rpt/drfine+1
                    j = zpt/dzfine+1

					if( i > rmax-rmin .or. j > zmax-zmin ) then
						write(*,*) "ERROR! TRYING to access:"
						write(*,*) "i,j", i, j
						write(*,*) "max i:", rmax-rmin
						write(*,*) "max j:", zmax-zmin
						write(*,*) rpt, "rpt"
						write(*,*) zpt, "zpt"
						write(*,*) drfine, "drfine"
						write(*,*) dzfine, "dzfine"
						write(*,*) rmax, "rmax"
						write(*,*) rmin, "rmin"
						write(*,*) zmax, "zmax"
						write(*,*) zmin, "zmin"
						stop 777
					end if

					q = findcell(i,j)%proc_i
					r = findcell(i,j)%proc_j 
					s = findcell(i,j)%cl_index 

                    processors(q,r)%collcells(s)%mask = maskcode

                end if
                ! bottom of in-bounds if

            end do
        end do
        ! end of circle drawing


        ! now draw the line halfway to the next point

        nline = .5*d/drho+1
        ds = .5*d/real(nline)

        do ii=1,nline ! step along the line segment
        
            do jj=-numsteps,numsteps ! paint on both sides of the line
                zpt=z0+ii*ds*wz+jj*drho*gz
                rpt=r0+ii*ds*wx+jj*drho*gx

                ! top of in-bounds if
                if(rpt.lt.(rmax*drfine).and.rpt.gt.(rmin*drfine).and. &
                    zpt.gt.(zmin*dzfine).and.zpt.lt.(zmax*dzfine)) then

                    i = rpt/drfine+1
                    j = zpt/dzfine+1

					if( i > rmax-rmin .or. j > zmax-zmin ) then
						write(*,*) "ERROR! TRYING to access:"
						write(*,*) "i,j", i, j
						write(*,*) "max i:", rmax-rmin
						write(*,*) "max j:", zmax-zmin
						write(*,*) rpt, "rpt"
						write(*,*) zpt, "zpt"
						write(*,*) drfine, "drfine"
						write(*,*) dzfine, "dzfine"
						write(*,*) rmax, "rmax"
						write(*,*) rmin, "rmin"
						write(*,*) zmax, "zmax"
						write(*,*) zmin, "zmin"
						stop 888
					end if

					q = findcell(i,j)%proc_i
					r = findcell(i,j)%proc_j 
					s = findcell(i,j)%cl_index 

                    processors(q,r)%collcells(s)%mask = maskcode

                end if
                ! bottom of in-bounds if

            end do
        end do

    end subroutine painter



    subroutine mask_metal(polygon)
!This subroutine uses a "crossing number" algorithm to determine whether or not a
!collision cell is inside a polygon (inside metal)
!The collision cell is then given a negative mask code that corresponds to the
!polygon number that it is inside
!mask = -1 is inside polygon number 1
!mask = -2 is inside polygon number 2, and so on
        use data_mod
        use types_mod
        implicit none
       
         
        ! Build a type for storing information about line segments
        type lineseg_type
           integer, allocatable, dimension(:) :: flag, pointflag, skipflag, rangeflag
           real(8), allocatable, dimension(:) :: rmin, rmax, zmin,zmax, rstart
           real(8), allocatable, dimension(:) :: slope, intercept 
           integer Nlines

        end type lineseg_type
         
        type(polygon_type) polygon(npoly)
        type(lineseg_type) lines(npoly)

        real(8) rr,zz, zcross
        integer i,j,k,m,ii, crosses, maskint
        integer ibefore,iafter,nextswitch

       !First, load the needed linesegment information
       do m=1,npoly

        !First, find the number of lines in the polygon
        if(polygon(m)%polytype.eq.1) lines(m)%Nlines=polygon(m)%numpoints
        if(polygon(m)%polytype.eq.0) lines(m)%Nlines=polygon(m)%numpoints 
        !open polygons haven't been tested - may cause errors

           !Now allocate the variables that depend on the number of linesegments.
           allocate(lines(m)%rmin(lines(m)%Nlines),lines(m)%rmax(lines(m)%Nlines),lines(m)%flag(lines(m)%Nlines),lines(m)%zmax(lines(m)%Nlines))
           allocate(lines(m)%zmin(lines(m)%Nlines),lines(m)%slope(lines(m)%Nlines),lines(m)%intercept(lines(m)%Nlines))
           allocate(lines(m)%skipflag(lines(m)%Nlines),lines(m)%pointflag(lines(m)%Nlines),lines(m)%rstart(lines(m)%Nlines))
           allocate(lines(m)%rangeflag(lines(m)%Nlines))
           !Now run over each lines segment to convert it to slope-intercept form
           do i=1,lines(m)%Nlines-1
                lines(m)%slope(i)=(polygon(m)%points(i+1)%r-polygon(m)%points(i)%r)*drfine/((polygon(m)%points(i+1)%z-polygon(m)%points(i)%z)*dzfine)
                lines(m)%intercept(i)=polygon(m)%points(i)%r*drfine-lines(m)%slope(i)*polygon(m)%points(i)%z*dzfine

                !Verticle and horizontal lines are given special cases, so flag
                !them
                lines(m)%flag(i)=0 !Set flag generally to zero, the correct for exceptions
                if(polygon(m)%points(i+1)%z.eq.polygon(m)%points(i)%z) lines(m)%flag(i)=1 !Verticle line
                if(polygon(m)%points(i+1)%r.eq.polygon(m)%points(i)%r) lines(m)%flag(i)=2 !Horizontal line

                !Now find the endpoints sorted by size, not input order
                lines(m)%rmax(i)=max(polygon(m)%points(i)%r*drfine,polygon(m)%points(i+1)%r*drfine)             
                lines(m)%rmin(i)=min(polygon(m)%points(i)%r*drfine,polygon(m)%points(i+1)%r*drfine)             
                lines(m)%zmin(i)=min(polygon(m)%points(i)%z*dzfine,polygon(m)%points(i+1)%z*dzfine)
                lines(m)%zmax(i)=max(polygon(m)%points(i)%z*dzfine,polygon(m)%points(i+1)%z*dzfine)
                !Find the starting radius of the line. It is used later to check
                !if the line starts at a radius that falls inside the collision
                !cell
                lines(m)%rstart(i)=dble(polygon(m)%points(i)%r*drfine)

            end do !End of loop over line segments
        !Now grab the last line segment from the last point to the first            
           
                lines(m)%slope(lines(m)%Nlines)=(polygon(m)%points(1)%r*drfine-polygon(m)%points(lines(m)%Nlines)%r*drfine)/(polygon(m)%points(1)%z*dzfine-polygon(m)%points(lines(m)%Nlines)%z*dzfine)
                lines(m)%intercept(lines(m)%Nlines)=polygon(m)%points(lines(m)%Nlines)%r*drfine-lines(m)%slope(lines(m)%Nlines)*polygon(m)%points(lines(m)%Nlines)%z*dzfine

                !Verticle and horizontal lines are given special cases, so flag
                !them
                lines(m)%flag(lines(m)%Nlines)=0 !Set flag generally to zero, then correct for exceptions
                if(polygon(m)%points(1)%z.eq.polygon(m)%points(lines(m)%Nlines)%z) lines(m)%flag(lines(m)%Nlines)=1 !Verticle line
                if(polygon(m)%points(1)%r.eq.polygon(m)%points(lines(m)%Nlines)%r) lines(m)%flag(lines(m)%Nlines)=2 !Horizontal line

                !Now find the endpoints sorted by size, not input order
                lines(m)%rmax(lines(m)%Nlines)=max(polygon(m)%points(1)%r*drfine,polygon(m)%points(lines(m)%Nlines)%r*drfine)             
                lines(m)%rmin(lines(m)%Nlines)=min(polygon(m)%points(1)%r*drfine,polygon(m)%points(lines(m)%Nlines)%r*drfine)             
                lines(m)%zmin(lines(m)%Nlines)=min(polygon(m)%points(1)%z*dzfine,polygon(m)%points(lines(m)%Nlines)%z*dzfine)
                lines(m)%zmax(lines(m)%Nlines)=max(polygon(m)%points(1)%z*dzfine,polygon(m)%points(lines(m)%Nlines)%z*dzfine)
                !Find the starting radius of the line. It is used later to check
                !if the line starts at a radius that falls inside the collision
                !cell
                lines(m)%rstart(lines(m)%Nlines)=dble(polygon(m)%points(lines(m)%Nlines)%r*drfine)

             !This section assigns a flag to each endpoint so that the metal
             !detection algorithm will know what to do with the point if it falls
             !inside a collision cell.
             !The assignment is made by looking at the points before and after
             !the current point. It the before and after points are both on the
             !same side of the current point, then none of the line segments
             !will cross through the collision cell and the point can be
             !ignored.  If the before and after points are on opposite sides of
             !the center point, then the line segments connecting them will pass
             !through the collision cell, and should be counted as an
             !intersection
             !For horizontal lines, the algorithm looks for the points before
             !and after the line to determine what type of point it is, and the
             !information is flagged on the last point of the line.

             do i=1,lines(m)%Nlines
                  nextswitch=0 !Used to search along horizontal lines. 
                  ibefore=i-1 !Previous point on the polynomial
                  !If ibefore reaches past the first point, have it go to the
                  !other end of the array
                  if(ibefore.lt.1) ibefore=lines(m)%Nlines+ibefore
                  iafter=i+1 !Next point
                  !if reaching past the end of the array, go back to the first
                  !point 
                  if(iafter.gt.lines(m)%Nlines) iafter=iafter-lines(m)%Nlines
                  

                  !Now check for horizontal segments.  Since they run parrallel
                  !to the line crossing ray, we only care about what happens at
                  !the end point, and can ignore any other points on the line
                  if(polygon(m)%points(i)%r.eq.polygon(m)%points(iafter)%r) then
                        lines(m)%pointflag(i)=0
                        cycle 
                  !Tell the last point on the line what the points before and
                  !after the line look like.
                  !If we've gotten to here, the iafter point will not be on the
                  !same line
                  else if (polygon(m)%points(i)%r.eq.polygon(m)%points(ibefore)%r) then 
                        ! Now find the point just before the line started
                        do while (nextswitch.eq.0)
                            !reach back one more point    
                            ibefore=ibefore-1
                            if(ibefore.lt.1) ibefore=lines(m)%Nlines+ibefore
                            !Now see if we've gone back far enough to get off
                            !the horizontal line.  If we have, flip the search
                            !switch
                            if(polygon(m)%points(i)%r.ne.polygon(m)%points(ibefore)%r) nextswitch=1
                        end do
                  end if

                  !If we reach this point, we are either on a normal point (not
                  !in a horziontal line) or at the endpoint of a horizontal line
                  !and ibefore is set to the point before the line started.
                  
                  !If before and after points are on the same side:
                  if(sign(1,nint(polygon(m)%points(i)%r-polygon(m)%points(iafter)%r)).eq. &
                     sign(1,nint(polygon(m)%points(i)%r-polygon(m)%points(ibefore)%r))) &
                     lines(m)%pointflag(i)=0
                  !if the before and after points are on opposite sides:
                  if(sign(1,nint(polygon(m)%points(i)%r-polygon(m)%points(iafter)%r)).ne. &
                     sign(1,nint(polygon(m)%points(i)%r-polygon(m)%points(ibefore)%r))) &
                     lines(m)%pointflag(i)=1

             end do
                   
 
       end do ! End of loop over polynomials 


       !Now run over every cell and count how many times a horizontal ray that
       !starts in the center of the cell and points in the positive z direction
       !crosses the line segments of each polynomial.  If the ray crosses an
       !even number of times, the point is outside the polynomial. If it crosses
       !an odd number of times, it is inside.

       !Now loop over every cell:

        do i = 1,num_proc_r
            do j = 1,num_proc_z
                do k = 1,processors(i,j)%nc
              
                    !First, find the center of the cell
                    rr=dble(processors(i,j)%collcells(k)%rmin+processors(i,j)%collcells(k)%rmax)/2d0*drfine
                    zz=dble(processors(i,j)%collcells(k)%zmin+processors(i,j)%collcells(k)%zmax)/2d0*dzfine

                   !Now flag line segments that have endpoints between rmin and
                   !rmax so the algorithm knows to skip them
                   !First, zero the old flag
                   do m=1,npoly
                   do ii=1,lines(m)%Nlines
                         lines(m)%skipflag(ii)=0
                   end do
                   end do
                   !Now check for endpoints
                   do m=1,npoly
                      do ii=1,lines(m)%Nlines
                         !Check if the starting point is between rmin and rmax
                         ! That mean the ray intersection aglorithm could have
                         ! problems from the horizontal ray being too close to
                         ! the end point.  Each point has already been flagged
                         ! with information on what to do in this case, so just
                         ! load a variable that tells the algortihm to skip the
                         ! line before and after the endpoint that falls between
                         ! rmin and rmax
                         if(lines(m)%rstart(ii).gt.dble(processors(i,j)%collcells(k)%rmin)*drfine.and.&
                            lines(m)%rstart(ii).lt.dble(processors(i,j)%collcells(k)%rmax)*drfine) then
                            !Setting this to 1 lets the algorithm know that the
                            !endpoint falls between rmin and rmax.  rangeflag=0
                            !tells the algorithm that the point does not fall in
                            !this range
                            lines(m)%rangeflag(ii)=1
                            !This flag tells the algorithm to not search for
                            !intersections on the lines segments next to a point
                            !whose r value falls between rmin and rmax.  The
                            !endpoint that is being checked has the same index
                            !as the lines segment that follows it, so the
                            !algorithm needs to know to skip the lines segment
                            !with the same index as the point, and the line
                            !segment before it.
                            
                            lines(m)%skipflag(ii)=1
                            !Now reach back to the line before. If this reaches
                            !past zero, loop back to the end of the array.
                            !Nothing special has to be done if we are moving
                            !along a horizontal line. Reaching back along a line
                            !will just set the variable to one again.
                            ibefore=ii-1
                            if(ibefore.lt.1) ibefore=lines(m)%Nlines-ibefore
                            lines(m)%skipflag(ibefore)=1
                           
                         !Now make sure that there are no more rangeflags set
                         !from the previous cell
                         else
                            lines(m)%rangeflag(ii)=0
                         end if
                      end do !end of loop over lines segments
                   end do !end of loop over polynomials
 

                   !Output rangeflag and skipflag for debugging                   
!                  if(rr.gt.1.30684e-03.and.rr.lt.1.308918e-03) then
!                        write(*,*) '##',rr,zz,dble(processors(i,j)%collcells(k)%rmin),dble(processors(i,j)%collcells(k)%rmax)
!                  do m=1,npoly
!                     do ii=1,lines(m)%Nlines
!                        write(*,*) '     ',lines(m)%rstart(ii),nint(lines(m)%rstart(ii)/drfine/20),lines(m)%rangeflag(ii),lines(m)%skipflag(ii)
!                     end do
!                  end do
!                  stop 666

!                  end if 
                      
                    

                    !Now run over each endpoint and check for crossings
                    do m=1,npoly
                       maskint=-m  !This sets the mask for each polynomial
                       crosses=0  !reset counter for each polynomial
                       !Now go over each line segment of the polynomial and
                       !check for intersections
                       
                       
                       do ii=1,lines(m)%Nlines
                          !First, check for special points
                          if(lines(m)%skipflag(ii).eq.1) then
                             !Count crosses from special points
                             !Shape patterns that should have an intersection
                             !are given pointflag=1, shape patterns that
                             !shouldn't are given pointflag=0
                             if(lines(m)%rangeflag(ii).eq.1) crosses=crosses+lines(m)%pointflag(ii) 
                             !Skip the intersection algorithm
                             cycle
                          end if

                          !Now check for intersections on remaining line
                          !segments
                          !Check if the cell is in the right range of r
                          if(rr.gt.lines(m)%rmin(ii).and.rr.lt.lines(m)%rmax(ii)) then

                               !If the segment is horizontal, they cannot intersect, so only worry about flag=0,1 
                               !If the line is verticle, greater in z, they will intersect
                               if(lines(m)%flag(ii).eq.1.and.zz.lt.lines(m)%zmin(ii)) crosses=crosses+1
                               !For the other types of lines:
                               if(lines(m)%flag(ii).eq.0) then
                                  !Find the z value of the intersection
                                  zcross=(rr-lines(m)%intercept(ii))/lines(m)%slope(ii)
                                  
                                  !Now check if the intersection is greater in z
                                  !than then center of the cell
                                  
                                  if(zcross.gt.zz) crosses=crosses+1
            !                     end if
                               end if

                          end if

                       end do !End of loop over line segments
                       !At this point, all of the intersections between the ray
                       !and the current polynomial have been counted.
                       !Now check and see if the number of crossings is odd, and
                       !load the mask accordingly
                       if(mod(crosses,2).eq.1) processors(i,j)%collcells(k)%mask=maskint

                    end do !End of loop over polynomials

                end do !End of loop running over all cells
            end do !End of loop runing over all cells
        end do !End of loop running over all cells


    end subroutine mask_metal


! THIS FUNCTION IS UNUSED AND UNFINISHED:

!	! this function determines the proper mask code for a particular collision cell, 
!	! given the cell's position, and given the geometry.
!	integer function cell_mask( bottom, top, left, right, drfine, dzfine, npoly, polygons )
!		use types_mod
!		implicit none

!		real(8) bottom, top, left, right  ! in units of drfine/dzfine
!		real(8) drfine, dzfine
!		integer npoly
!		type(polygon_type) polygons


!		cell_mask = 0


!	end function cell_mask



end module mask_calc_mod
