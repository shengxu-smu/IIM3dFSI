c-----------------------------------------------------------------------
c
      subroutine lagrange_interextrapolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer isphr,nijk,nn,id,jd,kd,ic,ie,jc,je,kc,ke
      integer iu,ju,ku,iv,jv,kv,iw,jw,kw
      parameter(nijk=2)
      real*8 sn,xx,yy,zz,xn,yn,zn,foo,gacobi
      real*8 uu(2),vv(2),ww(2)
      real*8 hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)

      sn=1.01d0*sqrt(dx*dx+dy*dy+dz*dz)

      DO m=1,ms

      isphr=i01(m)

      do m2=0,ns2-1
        do m1=0,ns1

          gacobi=fst(m1,m2,10,m)

          do n=-1,1,2

            nn=(n+3)/2
            xn=dble(n)*fst(m1,m2,7,m)/gacobi
            yn=dble(n)*fst(m1,m2,8,m)/gacobi
            zn=dble(n)*fst(m1,m2,9,m)/gacobi
            id=int(sign(1.0,xn))
            jd=int(sign(1.0,yn))
            kd=int(sign(1.0,zn))
            xx=xs(m1,m2,m)+sn*xn
            yy=ys(m1,m2,m)+sn*yn
            zz=zs(m1,m2,m)+sn*zn

            i=int((xx-x0)/hdx)
            j=int((yy-y0)/hdy)
            k=int((zz-z0)/hdz)
            if(mod(i,2).eq.0) then
              ic=i/2+1
              ie=ic-1
            else
              ie=(i+1)/2
              ic=ie
            endif
            if(mod(j,2).eq.0) then
              jc=j/2+1
              je=jc-1
            else
              je=(j+1)/2
              jc=je
            endif
            if(mod(k,2).eq.0) then
              kc=k/2+1
              ke=kc-1
            else
              ke=(k+1)/2
              kc=ke
            endif

            iu=ie
            ju=jc
            ku=kc
            iv=ic
            jv=je
            kv=kc
            iw=ic
            jw=jc
            kw=ke
            if(id.lt.0.0d0) then
              iu=iu+1
              iv=iv+1
              iw=iw+1
            endif
            if(jd.lt.0.0d0) then
              ju=ju+1
              jv=jv+1
              jw=jw+1
            endif
            if(kd.lt.0.0d0) then
              ku=ku+1
              kv=kv+1
              kw=kw+1
            endif

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xe(iu+id*(i-1))
                  fi(i)=u(iu+id*(i-1),ju+jd*(j-1),ku+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                hj(j)=yc(ju+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              hk(k)=zc(ku+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,uu(nn),foo)

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xc(iv+id*(i-1))
                  fi(i)=v(iv+id*(i-1),jv+jd*(j-1),kv+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                hj(j)=ye(jv+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              hk(k)=zc(kv+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,vv(nn),foo)

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xc(iw+id*(i-1))
                  fi(i)=w(iw+id*(i-1),jw+jd*(j-1),kw+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                hj(j)=yc(jw+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              hk(k)=ze(kw+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,ww(nn),foo)

          enddo

          us(m1,m2,m)=0.5d0*(uu(1)+uu(2))-0.5d0*fr(m1,m2,1,m)*sn/gacobi
          vs(m1,m2,m)=0.5d0*(vv(1)+vv(2))-0.5d0*fr(m1,m2,2,m)*sn/gacobi
          ws(m1,m2,m)=0.5d0*(ww(1)+ww(2))-0.5d0*fr(m1,m2,3,m)*sn/gacobi

        enddo
      enddo

      do m1=0,ns1
        us(m1,ns2,m)=us(m1,0,m)
        vs(m1,ns2,m)=vs(m1,0,m)
        ws(m1,ns2,m)=ws(m1,0,m)
      enddo

      if(isphr.eq.0) then
        do m2=0,ns2
          us(ns1,m2,m)=us(0,m2,m)
          vs(ns1,m2,m)=vs(0,m2,m)
          ws(ns1,m2,m)=ws(0,m2,m)
        enddo
      endif

      ENDDO
       
      call filter(us,ms,8,8)
      call filter(vs,ms,8,8)
      call filter(ws,ms,8,8)

      return
      end


c-----------------------------------------------------------------------
