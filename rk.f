c-----------------------------------------------------------------------
c
      subroutine rk
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk
      real*8 fac

c  krk=1:
      krk=1
      fac=0.5d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      call old_save

      if(imove.eq.1) call surface_move(krk,fac)
      call pressure(fac)
      call rhs

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
          enddo
        enddo
      enddo 
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
          enddo
        enddo
      enddo

      call ubc
      call vbc
      call wbc

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            um(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            vm(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            wm(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo

c  krk=2:
      krk=2
      fac=0.5d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      if(imove.eq.1) call surface_move(krk,fac)
      call pressure(fac)
      call rhs

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
          enddo
        enddo
      enddo

      call ubc
      call vbc
      call wbc

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo

c  krk=3:
      krk=3
      fac=1.0d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      if(imove.eq.1) call surface_move(krk,fac)
      call pressure(fac)
      call rhs

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
          enddo
        enddo
      enddo

      call ubc
      call vbc
      call wbc

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo

c  krk=4:
      krk=4
      fac=1.0d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      if(imove.eq.1) call surface_move(krk,fac)
      call pressure(fac)
      call rhs

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            u(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            v(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
          enddo
        enddo
      enddo
      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            w(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo

      call ubc
      call vbc
      call wbc

      return
      end


c-----------------------------------------------------------------------
