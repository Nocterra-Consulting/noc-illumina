       subroutine twodin(nbx,nby,filename,bindata)
c read a real array from an unformatted binary file
c nbx,nby are inputs: the domain size taken from the topography file
c (see twodsize). The file header must have the same size.
       integer nbx,nby,i,j,nbxf,nbyf
       real bindata(nbx,nby)
       character*(*) filename
       open(unit=1,form='unformatted',file=filename,action='read')
         read(1) nbxf,nbyf
         if ((nbxf.ne.nbx).or.(nbyf.ne.nby)) then
          print*,'Domain size mismatch in file: ',filename
          print*,'File size is: ',nbxf,'x',nbyf
          print*,'Topography size is: ',nbx,'x',nby
          print*,'Computation aborted'
          stop 1
         endif
         do j=nby,1,-1
            do i=1,nbx
               read(1) bindata(i,j)
            enddo
         enddo
       close(unit=1)
       return
       end
c-----------------------------------------------------------------------
       subroutine twodsize(filename,nbx,nby)
c read only the nbx,nby header of an unformatted binary file
       integer nbx,nby
       character*(*) filename
       open(unit=1,form='unformatted',file=filename,action='read')
         read(1) nbx,nby
       close(unit=1)
       if ((nbx.lt.1).or.(nby.lt.1)) then
         print*,'Invalid domain size in file: ',filename
         print*,'Size read: ',nbx,'x',nby
         print*,'Computation aborted'
         stop 1
       endif
       return
       end
