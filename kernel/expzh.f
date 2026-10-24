c-----------------------------------------------------------------------
c
c=======================================================================
c Function expzh
c
c Returns dexp(-1.*z/h), the relative density at the altitude z for the
c scale height h, from a cache of the last four altitudes per profile.
c The 2nd order scattering loop asks for the same few altitudes (the
c source, the reflecting surface, the scattering voxel and the line of
c sight voxel) in every transmittance and scattering call; the cache
c returns the value computed by the same expression, so the results are
c unchanged.
c
c is = 1 molecules, 2 aerosols, 3 layer
c-----------------------------------------------------------------------
      real*8 function expzh(z,h,is)
      real*8 z,h
      integer is,i
      real*8 kz(4,3),kh(4,3),ez(4,3)
      integer nx(3)
      save kz,kh,ez,nx
!$omp threadprivate(kz,kh,ez,nx)
      data kz /12*-1.d30/, kh /12*0.d0/, nx /3*0/
      do i=1,4
        if ((kz(i,is).eq.z).and.(kh(i,is).eq.h)) then
          expzh=ez(i,is)
          return
        endif
      enddo
      nx(is)=mod(nx(is),4)+1
      kz(nx(is),is)=z
      kh(nx(is),is)=h
      ez(nx(is),is)=dexp(-1.*z/h)
      expzh=ez(nx(is),is)
      return
      end
c-----------------------------------------------------------------------
c
c=======================================================================
c Function cosz
c
c Returns dcos(a) for the last angle asked. transmitm, transmita and
c transmitl are called in turn with the same zenith angle, so two of
c the three cosines come from this cache. The value is unchanged.
c-----------------------------------------------------------------------
      real*8 function cosz(a)
      real*8 a
      real*8 ka,ca
      save ka,ca
!$omp threadprivate(ka,ca)
      data ka /-1.d30/
      if (a.ne.ka) then
        ca=dcos(a)
        ka=a
      endif
      cosz=ca
      return
      end
c-----------------------------------------------------------------------
c
c=======================================================================
c Function scatf
c
c Returns 1.-dexp(lg*dexp(-1.*z/h)*un/h), the scattering factor of
c diffusion, from a cache of the last four calls per profile. The
c factor depends only on the altitude and on run constants; the
c scattering calls of one scattering voxel share the altitude, so most
c factors come from the cache. The value is unchanged.
c
c is = 1 molecules, 2 aerosols, 3 layer
c-----------------------------------------------------------------------
      real*8 function scatf(z,lg,un,h,is)
      real*8 z,lg,un,h
      integer is,i
      real*8 kz(4,3),kl(4,3),ku(4,3),kh(4,3),fz(4,3)
      real*8 expzh
      integer nx(3)
      save kz,kl,ku,kh,fz,nx
!$omp threadprivate(kz,kl,ku,kh,fz,nx)
      data kz /12*-1.d30/, kl /12*0.d0/, ku /12*0.d0/, kh /12*0.d0/
      data nx /3*0/
      do i=1,4
        if ((kz(i,is).eq.z).and.(kl(i,is).eq.lg).and.
     +  (ku(i,is).eq.un).and.(kh(i,is).eq.h)) then
          scatf=fz(i,is)
          return
        endif
      enddo
      nx(is)=mod(nx(is),4)+1
      kz(nx(is),is)=z
      kl(nx(is),is)=lg
      ku(nx(is),is)=un
      kh(nx(is),is)=h
      fz(nx(is),is)=1.-dexp(lg*expzh(z,h,is)*un/h)
      scatf=fz(nx(is),is)
      return
      end
