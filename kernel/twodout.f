       subroutine twodout(nbx,nby,filename,bindata)
c write a real array to an unformatted binary file
       integer nbx,nby,i,j
       real*8 bindata(nbx,nby)
       character*(*) filename
       open(unit=1,form='unformatted',file=filename,action='write')
         write(1) nbx,nby
         do j=nby,1,-1
            do i=1,nbx
               write(1) SNGL(bindata(i,j))
            enddo
         enddo
       close(unit=1)
       return
       end
