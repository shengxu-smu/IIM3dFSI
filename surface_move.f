c-----------------------------------------------------------------------
c
      subroutine surface_move(krk,fac)
      include 'parameter.inc'
      include 'surface.inc'
      integer krk,iprescribe
      real*8 fac

      iprescribe=1

c     call lagrange_interpolate
c     call lagrange_extrapolate
      call lagrange_interextrapolate

      DO m=1,ms

      do m2=0,ns2
        do m1=0,ns1
          usrk(krk,m1,m2,m)=us(m1,m2,m)
          vsrk(krk,m1,m2,m)=vs(m1,m2,m)
          wsrk(krk,m1,m2,m)=ws(m1,m2,m)
          if(krk.eq.4) then
            xs(m1,m2,m)=xsn(m1,m2,m)+(dt/6.0d0)*(usrk(1,m1,m2,m)+
     .        2.0d0*(usrk(2,m1,m2,m)+usrk(3,m1,m2,m))+usrk(4,m1,m2,m))
            ys(m1,m2,m)=ysn(m1,m2,m)+(dt/6.0d0)*(vsrk(1,m1,m2,m)+
     .        2.0d0*(vsrk(2,m1,m2,m)+vsrk(3,m1,m2,m))+vsrk(4,m1,m2,m))
            zs(m1,m2,m)=zsn(m1,m2,m)+(dt/6.0d0)*(wsrk(1,m1,m2,m)+
     .        2.0d0*(wsrk(2,m1,m2,m)+wsrk(3,m1,m2,m))+wsrk(4,m1,m2,m))
          else
            xs(m1,m2,m)=xsn(m1,m2,m)+fac*dt*us(m1,m2,m)
            ys(m1,m2,m)=ysn(m1,m2,m)+fac*dt*vs(m1,m2,m)
            zs(m1,m2,m)=zsn(m1,m2,m)+fac*dt*ws(m1,m2,m)
          endif
        enddo
      enddo

      ENDDO

      if(iprescribe.eq.1) call surface_prescribe(krk,fac)

      return
      end


c-----------------------------------------------------------------------
