c-----------------------------------------------------------------------
c
c=======================================================================
c Routine transmitl 
c
c Determine la transmittance des aerosols sur un parcours entre les 
c cellules (x_i,y_i,z_i) et (x_f,y_f,z_f)
c Est fonction de l'epaisseur optique des aerosols
c Recoit des longueurs d'ondes en nanometre et les transforme en microns.
c La pression doit etre en KPa
c Retourne la transmittance transl
c
c
c pour utilisation avec Illumina 
c-----------------------------------------------------------------------
c   
c    Copyright (C) 2021  Martin Aube
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c    Contact: martin.aube@cegepsherbrooke.qc.ca
c
c
      subroutine transmitl(angz,z_i,z_f,distd,hlay,transl,tranal)
      real*8 transl                                                         ! Declaration des variables.
      real*8 tranal                                                         ! vertical transmittance of the complete atmosphere (aerosols)
      real*8 angz,hlay
      real*8 z_i,z_f,z1,z2
      real*8 expzh,cosz
      real*8 distd
      real*8 lgtal,k_lgtal
      save lgtal,k_lgtal
!$omp threadprivate(lgtal,k_lgtal)
      data k_lgtal /-1.d0/


      if (z_i.gt.z_f) then
        z2=z_i
        z1=z_f
      else
        z1=z_i
        z2=z_f
      endif
c lgtal = dlog(tranal), cached: the argument is constant for a run, and
c dlog of the same argument gives the same value (results unchanged)
      if (tranal.ne.k_lgtal) then
        lgtal=dlog(tranal)
        k_lgtal=tranal
      endif
      if (z1.ne.z2) then    
        transl=dexp((lgtal/dabs(cosz(angz)))*(expzh(z1,hlay,3)-
     +  expzh(z2,hlay,3)))
      else
        transl=dexp((lgtal)*expzh(z1,hlay,3)*distd/hlay)  
      endif
         if (transl.eq.0.) then
            print*,'ERREUR transl - no transmission',z_i,z_f,angz
         endif
         if (transl.gt.1.) then 
            print*,'ERREUR avec transa',transl,z_i,z_f,angz
            stop 1
         endif
      return
      end
