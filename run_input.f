c----------------------------------------------------------------------
c
      subroutine run_input
      include 'parameter.inc'
      include 'surface.inc'

      open(unit=8,file='DAT/input.run',status='old')
      rewind 8
      read(8,*)
      read(8,*)nstep,icfl,iread,iwrite,itout,iplot,ianimation
      read(8,*)
      read(8,*)
      read(8,*)cflc,cflv,dtcfl,dtfix,tout
      close(8)

      open(unit=88,file='DAT/object.run',status='old')
      rewind 88
      read(88,*)
      read(88,*)
      read(88,*)
      do l=1,ms
        read(88,*)
        read(88,*)
        read(88,*)sk0(l),sk1(l),sk2(l),sk3(l)
        read(88,*)i01(l),ea(l),eb(l),ec(l),
     .            xsc0(l),ysc0(l),zsc0(l),phase(l)
        read(88,*)pitch0(l),yaw0(l),roll0(l),
     .            phasep(l),phasey(l),phaser(l)
      enddo
      close(88)

      open(unit=27,file='DAT/cxcycz.dat',status='unknown')

      return
      end


c-----------------------------------------------------------------------
