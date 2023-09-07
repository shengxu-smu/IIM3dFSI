c-----------------------------------------------------------------------
c
      subroutine surface_parametrization
      include 'parameter.inc'
      include 'surface.inc'
      real*8 angle1,angle2,rr,rx,ry,rz
      real*8 c,s,e,rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz

      open(unit=21,file='DAT/alfa10.dat',status='unknown')
      open(unit=22,file='DAT/alfa11.dat',status='unknown')
      open(unit=23,file='DAT/alfa20.dat',status='unknown')
      do m1=0,ns1
        alfa10(m1)=dalfa10*dble(m1)
        alfa11(m1)=dalfa11*dble(m1)+0.5d0*dalfa11
c       alfa11(m1)=dacos(1.0d0-2.0d0*alfa11(m1)/pi)
        write(21,100)alfa10(m1)
        write(22,100)alfa11(m1)
      enddo
      do m2=0,ns2
        alfa20(m2)=dalfa20*dble(m2)
        write(23,100)alfa20(m2)
      enddo
      close(21)
      close(22)
      close(23) 

      open(unit=24,file='DAT/xs0.dat',status='unknown')
      open(unit=25,file='DAT/ys0.dat',status='unknown')
      open(unit=26,file='DAT/zs0.dat',status='unknown')

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

      rxx=e*rx*rx+c
      rxy=e*rx*ry-s*rz
      rxz=e*rx*rz+s*ry
      ryx=e*rx*ry+s*rz
      ryy=e*ry*ry+c
      ryz=e*ry*rz-s*rx
      rzx=e*rx*rz-s*ry
      rzy=e*ry*rz+s*rx
      rzz=e*rz*rz+c

      if(i01(m).eq.1) then
c       rr=(ea(m)*eb(m)*ec(m))**(1.0d0/3.0d0)
        do m2=0,ns2
          angle2=alfa20(m2)
          do m1=0,ns1        
            angle1=alfa11(m1)
            xss(m1,m2,m)=ea(m)*dsin(angle1)*dcos(angle2)
            yss(m1,m2,m)=eb(m)*dsin(angle1)*dsin(angle2)
            zss(m1,m2,m)=ec(m)*dcos(angle1)
c           xss(m1,m2,m)=ea(m)*(2.0d0*dsqrt(angle1*(pi-angle1))/pi)*
c     .                  dcos(angle2)
c           yss(m1,m2,m)=eb(m)*(2.0d0*dsqrt(angle1*(pi-angle1))/pi)*
c     .                  dsin(angle2)
c           zss(m1,m2,m)=ec(m)*(2.0d0/pi)*(pi/2.0d0-angle1)
            xse(m1,m2,m)=xsc(m)+rxx*xss(m1,m2,m)+
     .                          rxy*yss(m1,m2,m)+rxz*zss(m1,m2,m)
            yse(m1,m2,m)=ysc(m)+ryx*xss(m1,m2,m)+
     .                          ryy*yss(m1,m2,m)+ryz*zss(m1,m2,m)
            zse(m1,m2,m)=zsc(m)+rzx*xss(m1,m2,m)+
     .                          rzy*yss(m1,m2,m)+rzz*zss(m1,m2,m)
c           xse(m1,m2,m)=xsc(m)+rr*dsin(angle1)*dcos(angle2)
c           yse(m1,m2,m)=ysc(m)+rr*dsin(angle1)*dsin(angle2)
c           zse(m1,m2,m)=zsc(m)+rr*dcos(angle1)
            xs(m1,m2,m)=xse(m1,m2,m)
            ys(m1,m2,m)=yse(m1,m2,m)
            zs(m1,m2,m)=zse(m1,m2,m)
          enddo
        enddo
      else
        do m2=0,ns2
          angle2=alfa20(m2)
          do m1=0,ns1
            angle1=alfa10(m1)
            xss(m1,m2,m)=ec(m)*dcos(angle1)
            yss(m1,m2,m)=(ec(m)*dsin(angle1)+ea(m))*dcos(angle2)
            zss(m1,m2,m)=(ec(m)*dsin(angle1)+eb(m))*dsin(angle2)
            xse(m1,m2,m)=xsc(m)+rxx*xss(m1,m2,m)+
     .                          rxy*yss(m1,m2,m)+rxz*zss(m1,m2,m)
            yse(m1,m2,m)=ysc(m)+ryx*xss(m1,m2,m)+
     .                          ryy*yss(m1,m2,m)+ryz*zss(m1,m2,m)
            zse(m1,m2,m)=zsc(m)+rzx*xss(m1,m2,m)+
     .                          rzy*yss(m1,m2,m)+rzz*zss(m1,m2,m)
            xs(m1,m2,m)=xse(m1,m2,m)
            ys(m1,m2,m)=yse(m1,m2,m)
            zs(m1,m2,m)=zse(m1,m2,m)
          enddo
        enddo
      endif

      do m1=0,ns1
        xs(m1,ns2,m)=xs(m1,0,m)
        ys(m1,ns2,m)=ys(m1,0,m)
        zs(m1,ns2,m)=zs(m1,0,m)
        xse(m1,ns2,m)=xse(m1,0,m)
        yse(m1,ns2,m)=yse(m1,0,m)
        zse(m1,ns2,m)=zse(m1,0,m)
        xss(m1,ns2,m)=xss(m1,0,m)
        yss(m1,ns2,m)=yss(m1,0,m)
        zss(m1,ns2,m)=zss(m1,0,m)
      enddo
      if(i01(m).eq.0) then
        do m2=0,ns2
          xs(ns1,m2,m)=xs(0,m2,m)
          ys(ns1,m2,m)=ys(0,m2,m)
          zs(ns1,m2,m)=zs(0,m2,m)
          xse(ns1,m2,m)=xse(0,m2,m)
          yse(ns1,m2,m)=yse(0,m2,m)
          zse(ns1,m2,m)=zse(0,m2,m)
          xss(ns1,m2,m)=xss(0,m2,m)
          yss(ns1,m2,m)=yss(0,m2,m)
          zss(ns1,m2,m)=zss(0,m2,m)
        enddo
      endif

      do m2=0,ns2
        write(24,100)(xs(m1,m2,m),m1=0,ns1)
        write(25,100)(ys(m1,m2,m),m1=0,ns1)
        write(26,100)(zs(m1,m2,m),m1=0,ns1)
      enddo

      ENDDO

      close(24)
      close(25)
      close(26)
100   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
