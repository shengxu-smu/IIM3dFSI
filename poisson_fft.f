c-----------------------------------------------------------------------
c
      subroutine poisson_fft
      include 'parameter.inc'
      include 'field.inc'
      include 'fft.inc'

      do k=1,nz

        do i=1,nx
          do j=1,ny
            ri(j,i)=o(i,j,k)*dx*dx
            if(k.eq.1) ri(j,i)=ri(j,i)+2.0d0*betaz*dz*pb(i,j)
            if(k.eq.nz) ri(j,i)=ri(j,i)-2.0d0*betaz*dz*pt(i,j)
          enddo
        enddo          

        do j=1,ny
          ri(j,1)=ri(j,1)+2.0d0*dx*pw(j,k)
          ri(j,nx)=ri(j,nx)-2.0d0*dx*pe(j,k)
        enddo
        
        do i=1,nx
          ri(1,i)=ri(1,i)+2.0d0*betay*dy*ps(i,k)
          ri(ny,i)=ri(ny,i)-2.0d0*betay*dy*pn(i,k)
        enddo

        call vcost(ny,nx,ri,wi,ny,wsavei)

        do j=1,ny
          do i=1,nx
            rj(i,j)=ri(j,i)
          enddo
        enddo

        call vcost(nx,ny,rj,wj,nx,wsavej)

        do j=1,ny
          do i=1,nx
            o(i,j,k)=rj(i,j)
          enddo
        enddo
  
      enddo

      do i=1,nx

        do k=1,nz
          do j=1,ny
            rk(j,k)=o(i,j,k)
          enddo
        enddo

        call vcost(ny,nz,rk,wk,ny,wsavek)

        do k=1,nz
          do j=1,ny
            if(i.eq.1.and.j.eq.1.and.k.eq.1) then
              rk(j,k)=0.0d0
            else
              rk(j,k)=0.5d0*rk(j,k)/(dcos(dble(i-1)*pi/dble(nx-1))+
     .                         betay*dcos(dble(j-1)*pi/dble(ny-1))+
     .                         betaz*dcos(dble(k-1)*pi/dble(nz-1))-
     .                         (1.0d0+betay+betaz))
            endif
          enddo
        enddo

        call vcost(ny,nz,rk,wk,ny,wsavek)

        do k=1,nz
          do j=1,ny
            o(i,j,k)=rk(j,k)
          enddo
        enddo

      enddo

      do k=1,nz

        do i=1,nx
          do j=1,ny
            ri(j,i)=o(i,j,k)
          enddo
        enddo

        call vcost(ny,nx,ri,wi,ny,wsavei)

        do j=1,ny
          do i=1,nx
            rj(i,j)=ri(j,i)
          enddo
        enddo

        call vcost(nx,ny,rj,wj,nx,wsavej)

        do j=1,ny
          do i=1,nx
            p(i,j,k)=rj(i,j)
          enddo
        enddo

      enddo

      return
      end


c-----------------------------------------------------------------------
