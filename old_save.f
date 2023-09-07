c-----------------------------------------------------------------------
c
      subroutine old_save
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

      if(imove.eq.1) then
        DO m=1,ms

         xscn(m)=xsc(m)
         yscn(m)=ysc(m)
         zscn(m)=zsc(m)
         xsctn(m)=xsct(m)
         ysctn(m)=ysct(m)
         zsctn(m)=zsct(m)

          do m2=0,ns2
            do m1=0,ns1
              xsn(m1,m2,m)=xs(m1,m2,m)
              ysn(m1,m2,m)=ys(m1,m2,m)
              zsn(m1,m2,m)=zs(m1,m2,m)
            enddo
          enddo
        ENDDO
      endif

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx
            un(i,j,k)=u(i,j,k)
            um(i,j,k)=u(i,j,k)
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=0,nx+1
            vn(i,j,k)=v(i,j,k)
            vm(i,j,k)=v(i,j,k)
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny+1
          do i=0,nx+1
            wn(i,j,k)=w(i,j,k)
            wm(i,j,k)=w(i,j,k)
          enddo
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          do i=1,nx
            dn(i,j,k)=data1(i,j,k)+data5(i,j,k)+data9(i,j,k)
          enddo
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------

