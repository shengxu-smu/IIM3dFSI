c-----------------------------------------------------------------------
c
      subroutine surface_prescribe(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      integer krk,motion
      real*8 fac,rx,ry,rz,rr,c,s,e,rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
      real*8 omega,tc,rxxt,rxyt,rxzt,ryxt,ryyt,ryzt,rzxt,rzyt,rzzt
      real*8 cz(ms),gbar,totm

      if(krk.eq.1) t0=t0+0.5d0*dt
      if(krk.eq.3) t0=t0+0.5d0*dt

      gbar=0.89331d0
      totm=(2.56d0-1.0d0)*(4.0d0/3.0d0)*pi*0.5d0*0.5d0*0.5d0

      do m=1,ms
        cz(m)=0.0d0
        do m2=0,ns2-1
          do m1=0,ns1-1
            cz(m)=cz(m)-0.25d0*(fz(m1,m2,m)*fst(m1,m2,10,m)+
     .                          fz(m1+1,m2,m)*fst(m1+1,m2,10,m)+
     .                          fz(m1,m2+1,m)*fst(m1,m2+1,10,m)+
     .                          fz(m1+1,m2+1,m)*fst(m1+1,m2+1,10,m))
           enddo
        enddo
        cz(m)=cz(m)*dalfa11*dalfa20
      enddo

      m=1

      zscrk(krk,m)=zsct(m)
      zsctrk(krk,m)=-gbar+cz(m)/totm

      if(krk.eq.4) then
        zsc(m)=zscn(m)+(dt/6.0d0)*(zscrk(1,m)+
     .          2.0d0*(zscrk(2,m)+zscrk(3,m))+zscrk(4,m))
        zsct(m)=zsctn(m)+(dt/6.0d0)*(zsctrk(1,m)+
     .          2.0d0*(zsctrk(2,m)+zsctrk(3,m))+zsctrk(4,m))
      else
        zsc(m)=zscn(m)+fac*dt*zscrk(krk,m)
        zsct(m)=zsctn(m)+fac*dt*zsctrk(krk,m)
      endif

c      omega=1.0d0
c      tc=1.0d0

c      rx=0.0d0
c      ry=-1.0d0
c      rz=0.0d0
c      rr=dsqrt(rx*rx+ry*ry+rz*rz)
c      rx=rx/rr
c      ry=ry/rr
c      rz=rz/rr

      DO m=1,ms

c      xsc(m)=xsc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3)
c      ysc(m)=ysc0(m)
c      zsc(m)=zsc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3)
c      pitch(m)=pitch0(m)
c      yaw(m)=yaw0(m)
c      roll(m)=roll0(m)+pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0)*
c     .                    (1.0d0-dexp(-t0/tc)))

c      xsct(m)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3)
c      ysct(m)=0.0d0
c      zsct(m)=1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3)
c      pitcht(m)=0.0d0
c      yawt(m)=0.0d0
c      rollt(m)=0.25d0*pi*(0.8d0*dcos(0.8d0*t0))*(1.0d0-dexp(-t0/tc))+
c     .         0.25d0*pi*(dsin(0.8d0*t0))*(dexp(-t0/tc)/tc)

c      xsctt(m)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3)
c      ysctt(m)=0.0d0
c      zsctt(m)=1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3)

      c=dcos(roll(m))
      s=dsin(roll(m))
      e=1.0d0-c

      rxx=e*rx*rx+c
      rxy=e*rx*ry-s*rz
      rxz=e*rx*rz+s*ry
      ryx=e*rx*ry+s*rz
      ryy=e*ry*ry+c
      ryz=e*ry*rz-s*rx
      rzx=e*rx*rz-s*ry
      rzy=e*ry*rz+s*rx
      rzz=e*rz*rz+c

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
          xse(m1,m2,m)=xsc(m)+rxx*xss(m1,m2,m)+
     .                        rxy*yss(m1,m2,m)+rxz*zss(m1,m2,m)
          yse(m1,m2,m)=ysc(m)+ryx*xss(m1,m2,m)+
     .                        ryy*yss(m1,m2,m)+ryz*zss(m1,m2,m)
          zse(m1,m2,m)=zsc(m)+rzx*xss(m1,m2,m)+
     .                        rzy*yss(m1,m2,m)+rzz*zss(m1,m2,m)
          use(m1,m2,m)=xsct(m)+rxxt*xss(m1,m2,m)+
     .                         rxyt*yss(m1,m2,m)+rxzt*zss(m1,m2,m)
          vse(m1,m2,m)=ysct(m)+ryxt*xss(m1,m2,m)+
     .                         ryyt*yss(m1,m2,m)+ryzt*zss(m1,m2,m)
          wse(m1,m2,m)=zsct(m)+rzxt*xss(m1,m2,m)+
     .                         rzyt*yss(m1,m2,m)+rzzt*zss(m1,m2,m)
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
