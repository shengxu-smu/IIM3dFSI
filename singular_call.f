c----------------------------------------------------------------------
c
      subroutine singular_call
      integer inc

      inc=1

      IF(inc.eq.1) THEN

      call surface_property
c     call surface_property_fd
      call singular_force

      call jc_rhs
c     call jc_rhs_fd
      call jc_first
      call jc_second

c     call derivative_check

      call euler_link
      call jc_velocity

c     call eulerlagrange_check

      call correction_interpolate
      call correction_strain

      call u_interpolate
      call u_strain
      call udu_surface

      call v_interpolate
      call v_strain
      call vdv_surface

      call w_interpolate
      call w_strain
      call wdw_surface        

      call jc_pressure
      call correction_difference

      ENDIF

      return
      end


c-----------------------------------------------------------------------
