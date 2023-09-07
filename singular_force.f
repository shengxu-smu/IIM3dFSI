c-----------------------------------------------------------------------
c
      subroutine singular_force
      include 'parameter.inc'
      include 'surface.inc'
      integer model(ms)
      real*8 xr1,yr1,zr1,xr2,yr2,zr2,xn,yn,zn,gacobi
      real*8 xr11,yr11,zr11,xr22,yr22,zr22,xr12,yr12,zr12
      real*8 sse,ssf,ssg,ssl,ssm,ssn,ssh
      real*8 foo(0:ns1,0:ns2)

      model(1)=2

      DO m=1,ms

      if(model(m).eq.1) then

      do m2=0,ns2
        do m1=0,ns1

          xr1=fst(m1,m2,1,m)
          yr1=fst(m1,m2,2,m)
          zr1=fst(m1,m2,3,m)
          xr2=fst(m1,m2,4,m)
          yr2=fst(m1,m2,5,m)
          zr2=fst(m1,m2,6,m)
          xn=fst(m1,m2,7,m)
          yn=fst(m1,m2,8,m)
          zn=fst(m1,m2,9,m)
          gacobi=fst(m1,m2,10,m)
          xr11=snd(m1,m2,1,m)
          yr11=snd(m1,m2,2,m)
          zr11=snd(m1,m2,3,m)
          xr22=snd(m1,m2,4,m)
          yr22=snd(m1,m2,5,m)
          zr22=snd(m1,m2,6,m)
          xr12=snd(m1,m2,7,m)
          yr12=snd(m1,m2,8,m)
          zr12=snd(m1,m2,9,m)
          sse=xr1*xr1+yr1*yr1+zr1*zr1
          ssf=xr1*xr2+yr1*yr2+zr1*zr2
          ssg=xr2*xr2+yr2*yr2+zr2*zr2
          ssl=xr11*xn+yr11*yn+zr11*zn
          ssm=xr12*xn+yr12*yn+zr12*zn
          ssn=xr22*xn+yr22*yn+zr22*zn
          ssh=(sse*ssn-2.0d0*ssf*ssm+ssg*ssl)/
     .        (2.0d0*(sse*ssg-ssf*ssf))/gacobi
          foo(m1,m2)=ssh
        enddo
      enddo

      call filter(foo,1,16,16)
      do m2=0,ns2
        do m1=0,ns1
          xn=fst(m1,m2,7,m)
          yn=fst(m1,m2,8,m)
          zn=fst(m1,m2,9,m)
          gacobi=fst(m1,m2,10,m)
          ssh=foo(m1,m2)
          fx(m1,m2,m)=sk0(m)*ssh*xn/gacobi
          fy(m1,m2,m)=sk0(m)*ssh*yn/gacobi
          fz(m1,m2,m)=sk0(m)*ssh*zn/gacobi
        enddo
      enddo

      endif

c

      if(model(m).eq.2) then
      do m2=0,ns2
        do m1=0,ns1

          fx(m1,m2,m)=sk1(m)*(use(m1,m2,m)-us(m1,m2,m))+
     .                sk2(m)*(xse(m1,m2,m)-xs(m1,m2,m))
          fy(m1,m2,m)=sk1(m)*(vse(m1,m2,m)-vs(m1,m2,m))+
     .                sk2(m)*(yse(m1,m2,m)-ys(m1,m2,m))
          fz(m1,m2,m)=sk1(m)*(wse(m1,m2,m)-ws(m1,m2,m))+
     .                sk2(m)*(zse(m1,m2,m)-zs(m1,m2,m))

        enddo
      enddo
      endif

      ENDDO
      
      return
      end


c-----------------------------------------------------------------------
