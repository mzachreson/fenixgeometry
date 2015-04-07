!------ onesource: combines all of the fenix10 source files into fenix10all.f90 --

   character*128 line,lineout
   character*30 filnam(20)
   integer i,j,k,m,n,ifiles


   ! define the source files to be written into fenix10all.f90

      filnam(3) = 'fgbuilder.f90'
      filnam(4) = 'fgdata.f90'
      filnam(2) = 'fgmain.f90'
      filnam(5) = 'fgmaskcalc.f90'
      filnam(1) = 'fgtypes.f90'
      filnam(6) = 'fgvolumecalc.f90'
      filnam(7) = 'fgwrite.f90'

      ifiles=7

   ! open fenix10all.f90
   open(unit=11,file='fenix10all.f90',status='unknown')

   ! now read each file one at a time, put a source file header into
   ! geometry10all.f90, and write the source file

   do i=1,ifiles

      open(unit=12,file=filnam(i),status='old')


      write(11,101) filnam(i)
101   format('Source file:   ',a30,'  *****************')
      write(11,*)

1     continue
      read(12,102,end=2) line
102   format(a128)
      k=len_trim(line)
      lineout=adjustl(line(1:k))
      write(11,102) lineout
      goto 1
2     continue

      close(unit=12)


   end do

   close(unit=11)

   stop
   end

