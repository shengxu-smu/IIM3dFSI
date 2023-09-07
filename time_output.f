c-----------------------------------------------------------------------
c
      subroutine time_output
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 tx,ty,bx,by,sz
      real*8 af(ms),cx(ms),cy(ms),cz(ms),vol(ms),ra(ms),rb(ms),rc(ms)

      do m=1,ms
        if(i01(m).eq.1) then
c         af(m)=pi*eb(m)*ec(m)
          af(m)=pi*ea(m)*eb(m)
        elseif(i01(m).eq.0) then
          af(m)=pi*(ea(m)+ec(m))*(eb(m)+ec(m))-
     .          pi*(ea(m)-ec(m))*(eb(m)-ec(m))
        endif
        cx(m)=0.0d0
        cy(m)=0.0d0
        cz(m)=0.0d0
        vol(m)=0.0d0
        do m2=0,ns2-1
          do m1=0,ns1-1
            cx(m)=cx(m)-0.25d0*(fx(m1,m2,m)*fst(m1,m2,10,m)+
     .                          fx(m1+1,m2,m)*fst(m1+1,m2,10,m)+
     .                          fx(m1,m2+1,m)*fst(m1,m2+1,10,m)+
     .                          fx(m1+1,m2+1,m)*fst(m1+1,m2+1,10,m))
            cy(m)=cy(m)-0.25d0*(fy(m1,m2,m)*fst(m1,m2,10,m)+
     .                          fy(m1+1,m2,m)*fst(m1+1,m2,10,m)+
     .                          fy(m1,m2+1,m)*fst(m1,m2+1,10,m)+
     .                          fy(m1+1,m2+1,m)*fst(m1+1,m2+1,10,m))
            cz(m)=cz(m)-0.25d0*(fz(m1,m2,m)*fst(m1,m2,10,m)+
     .                          fz(m1+1,m2,m)*fst(m1+1,m2,10,m)+
     .                          fz(m1,m2+1,m)*fst(m1,m2+1,10,m)+
     .                          fz(m1+1,m2+1,m)*fst(m1+1,m2+1,10,m))
            tx=xs(m1+1,m2,m)-xs(m1,m2,m)
            ty=ys(m1+1,m2,m)-ys(m1,m2,m)
            bx=xs(m1,m2+1,m)-xs(m1,m2,m)
            by=ys(m1,m2+1,m)-ys(m1,m2,m)
            sz=(zs(m1,m2,m)+zs(m1+1,m2,m)+zs(m1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
            tx=xs(m1+1,m2+1,m)-xs(m1,m2+1,m)
            ty=ys(m1+1,m2+1,m)-ys(m1,m2+1,m)
            bx=xs(m1+1,m2+1,m)-xs(m1+1,m2,m)
            by=ys(m1+1,m2+1,m)-ys(m1+1,m2,m)
            sz=(zs(m1+1,m2,m)+zs(m1+1,m2+1,m)+zs(m1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
           enddo
        enddo
        if(i01(m).eq.1) then 
          cx(m)=2.0d0*cx(m)*dalfa11*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa11*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa11*dalfa20/af(m)
          do m2=1,ns2-2
            tx=xs(0,m2+1,m)-xs(0,0,m)
            ty=ys(0,m2+1,m)-ys(0,0,m)
            bx=xs(0,m2+1,m)-xs(0,m2,m)
            by=ys(0,m2+1,m)-ys(0,m2,m)
            sz=(zs(0,0,m)+zs(0,m2,m)+zs(0,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)          
            tx=xs(ns1,0,m)-xs(ns1,m2+1,m)
            ty=ys(ns1,0,m)-ys(ns1,m2+1,m)
            bx=xs(ns1,m2+1,m)-xs(ns1,m2,m)
            by=ys(ns1,m2+1,m)-ys(ns1,m2,m)
            sz=(zs(ns1,0,m)+zs(ns1,m2,m)+zs(ns1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
          enddo
        else
          cx(m)=2.0d0*cx(m)*dalfa10*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa10*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa10*dalfa20/af(m)
        endif
        cx(m)=cx(m)+2.0d0*vol(m)*xsctt(m)/af(m)
        cy(m)=cy(m)+2.0d0*vol(m)*ysctt(m)/af(m)
        cz(m)=cz(m)+2.0d0*vol(m)*zsctt(m)/af(m)
 
        m1=int(ns1/2)
        m2=0
        ra(m)=dsqrt(xs(m1,m2,m)*xs(m1,m2,m)+ys(m1,m2,m)*ys(m1,m2,m)+
     .              zs(m1,m2,m)*zs(m1,m2,m))
        m1=int(ns1/2)
        m2=int(ns2/4)
        rb(m)=dsqrt(xs(m1,m2,m)*xs(m1,m2,m)+ys(m1,m2,m)*ys(m1,m2,m)+
     .              zs(m1,m2,m)*zs(m1,m2,m))
        m1=0
        m2=0
        rc(m)=dsqrt(xs(m1,m2,m)*xs(m1,m2,m)+ys(m1,m2,m)*ys(m1,m2,m)+
     .              zs(m1,m2,m)*zs(m1,m2,m))
      enddo
      write(27,100)t,(cx(m),cy(m),cz(m),vol(m),ra(m),rb(m),rc(m),
     .                zsc(m),zsct(m),m=1,ms)

 100  format(1x,200e16.6e4)

      return
      end


c-----------------------------------------------------------------------
