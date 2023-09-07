c-----------------------------------------------------------------------
c
      subroutine initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'fft.inc'
      real*8 tc,rx,ry,rz,rr,c,s,e
      real*8 rxxt,rxyt,rxzt,ryxt,ryyt,ryzt,rzxt,rzyt,rzzt

      call vrffti(ns1,wsave1)
      call vrffti(ns2,wsave2)
      call vrffti(2*ns1+2,wsave3)

      call vcosti(nx,wsavei)
      call vcosti(ny,wsavej)
      call vcosti(nz,wsavek)

      if(mod(ns2,2).ne.0) then
        write(*,*)'  !! erro: ns2 is not even !!'
        stop
      endif
      if(ms.gt.1) then
        write(*,*)'  !! note: keep objects away enough !!'
      endif

      t=0.0d0
      t0=0.0d0
      tc=1.0d0

      DO m=1,ms
         xsc(m)=xsc0(m)
         ysc(m)=ysc0(m)
         zsc(m)=zsc0(m)
c        xsc(m)=xsc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3)
c        ysc(m)=ysc0(m)
c        zsc(m)=zsc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3)
c        pitch(m)=pitch0(m)
c        yaw(m)=yaw0(m)
c        roll(m)=roll0(m)+pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0)*
c     .                      (1.0d0-dexp(-t0/tc)))
c
c        xsct(m)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3)
c        ysct(m)=0.0d0
c        zsct(m)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3)
c        pitcht(m)=0.0d0
c        yawt(m)=0.0d0
c        rollt(m)=0.25d0*pi*(0.8d0*dcos(0.8d0*t0))*(1.0d0-dexp(-t0/tc))+
c     .           0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc)
      ENDDO

      call surface_parametrization

      rx=0.0d0
      ry=-1.0d0
      rz=0.0d0
      rr=dsqrt(rx*rx+ry*ry+rz*rz)
      rx=rx/rr
      ry=ry/rr
      rz=rz/rr

      DO m=1,ms

      c=dcos(roll(m))
      s=dsin(roll(m))
      e=1.0d0-c

      rxxt=(s*rx*rx-s)*rollt(m)
      rxyt=(s*rx*ry-c*rz)*rollt(m)
      rxzt=(s*rx*rz+c*ry)*rollt(m)
      ryxt=(s*rx*ry+c*rz)*rollt(m)
      ryyt=(s*ry*ry-s)*rollt(m)
      ryzt=(s*ry*rz-c*rx)*rollt(m)
      rzxt=(s*rx*ry-c*ry)*rollt(m)
      rzyt=(s*ry*rz+c*rx)*rollt(m)
      rzzt=(s*rz*rz-s)*rollt(m)

      do m2=0,ns2
        do m1=0,ns1
          us(m1,m2,m)=0.0d0
          vs(m1,m2,m)=0.0d0
          ws(m1,m2,m)=0.0d0
          use(m1,m2,m)=xsct(m)+rxxt*xss(m1,m2,m)+
     .                         rxyt*yss(m1,m2,m)+rxzt*zss(m1,m2,m)
          vse(m1,m2,m)=ysct(m)+ryxt*xss(m1,m2,m)+
     .                         ryyt*yss(m1,m2,m)+ryzt*zss(m1,m2,m)
          wse(m1,m2,m)=zsct(m)+rzxt*xss(m1,m2,m)+
     .                         rzyt*yss(m1,m2,m)+rzzt*zss(m1,m2,m)
        enddo
      enddo
      ENDDO

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx
c           u(i,j,k)=1.0d0*dsin(yc(j))*dsin(zc(k))
            u(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=0,nx+1
c           v(i,j,k)=1.0d0*dsin(xc(i))*dsin(zc(k))
            v(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny+1
          do i=0,nx+1
c           w(i,j,k)=1.0d0*dsin(xc(i))*dsin(yc(j))
            w(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            p(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            d(i,j,k)=0.0d0
          enddo
        enddo
      enddo

c

      return
      end


c-----------------------------------------------------------------------
