       ! assume at the start that the segment between points 1 and 2 doesn't cross the cell defined by rbottom, rtop, zleft, and zright.
       ! The logic below will correct this if a crossing occurs

       does.. = .false.
       hitsbottom%r=-1.d0
       hitsbottom%z=-1.d0
       hitstop%r=-1.d0
       hitstop%z=-1.d0
       hitsright%r=-1.d0
       hitsright%z=-1.d0
       hitsleft%r=-1.d0
       hitsleft%z=-1.d0

       ! The calculation is based on the parametric form of the line between the two points segpt1%(r,z) and segpt2%(r,z)
       !  r(s) = (1-s)*r1+s*r2    z(s) = (1-s)*z1 + s*z2  with s between 0 and 1

       ! top of the cell

       if(segpt1%r.ne.segpt2%r) then ! for a top crossing r1=r2 means no crossing

          s=(rtop-segpt1%r)/(segpt2%r-segpt1%r)  ! find the s-value where the line crosses r=rtop
          z=(1.d0-s)*segpt1%z+s*segpt2%z ! find the z-value of the crossing

          ! check for a hit: s between 0 and 1 (on the segment) and z between zleft and zright (in the cell)

          if( s.gt.0.d0 .and s.lt.1.d0 . and. z.lt.zright .and. z.gt.zleft) then
             hitstop%r=rtop ! crossing at rtop
             hitstop%z=z ! z-value of the crossing at rtop
          end if

       end if

       ! the other sides of the cell are handled in a similar way

       ! bottom of the cell

       if(segpt1%r.ne.segpt2%r) then

          s=(rbottom-segpt1%r)/(segpt2%r-segpt1%r)
          z=(1.d0-s)*segpt1%z+s*segpt2%z

          ! check for a hit: s between 0 and 1 and z between zleft and zright
          if( s.gt.0.d0 .and s.lt.1.d0 . and. z.lt.zright .and. z.gt.zleft) then
             hitsbottom%r=rbottom
             hitsbottom%z=z
          end if

       end if

       ! right side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zright-segpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r

          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.0.d0 .and s.lt.1.d0 . and. r.lt.rtop .and. r.gt.rbottom) then
             hitsright%r=r
             hitsright%z=zright
          end if

       end if

       ! left side of the cell

       if(segpt1%z.ne.segpt2%z) then

          s=(zleft-segtpt1%z)/(segpt2%z-segpt1%z)
          r=(1.d0-s)*segpt1%r+s*segpt2%r

          ! check for a hit: s between 0 and 1 and r between rbottom and rtop
          if( s.gt.0.d0 .and s.lt.1.d0 . and. r.lt.rtop .and. r.gt.rbottom) then
             hitsleft%r=r
             hitsleft%z=zleft
          end if

       end if
