c-----------------------------------------------------------------------
c
      subroutine jc_velocity
      include 'parameter.inc'
      include 'surface.inc'
      real*8 signx,signy,signz

      DO m=1,ms

      do n=0,niejc(m)
        signx=sign(1.0,fiejc(10,n,m))
        signy=sign(1.0,fiejc(11,n,m))
        signz=sign(1.0,fiejc(12,n,m))
        fiejc(13,n,m)=signx*fiejc(13,n,m)
        fiejc(14,n,m)=signy*fiejc(14,n,m)
        fiejc(15,n,m)=signz*fiejc(15,n,m)
        fiejc(16,n,m)=signx*fiejc(16,n,m)
        fiejc(17,n,m)=signy*fiejc(17,n,m)
        fiejc(18,n,m)=signz*fiejc(18,n,m)
        fiejc(19,n,m)=signx*fiejc(19,n,m)
        fiejc(20,n,m)=signy*fiejc(20,n,m)
        fiejc(21,n,m)=signz*fiejc(21,n,m)
        fiejc(25,n,m)=signx*fiejc(25,n,m)
        fiejc(28,n,m)=signy*fiejc(28,n,m)
        fiejc(30,n,m)=signz*fiejc(30,n,m)
        fiejc(31,n,m)=signx*fiejc(31,n,m)
        fiejc(34,n,m)=signy*fiejc(34,n,m)
        fiejc(36,n,m)=signz*fiejc(36,n,m)
        fiejc(37,n,m)=signx*fiejc(37,n,m)
        fiejc(40,n,m)=signy*fiejc(40,n,m)
        fiejc(42,n,m)=signz*fiejc(42,n,m)
      enddo

      do n=0,nicje(m)
        signx=sign(1.0,ficje(10,n,m))
        signy=sign(1.0,ficje(11,n,m))
        signz=sign(1.0,ficje(12,n,m))
        ficje(13,n,m)=signx*ficje(13,n,m)
        ficje(14,n,m)=signy*ficje(14,n,m)
        ficje(15,n,m)=signz*ficje(15,n,m)
        ficje(16,n,m)=signx*ficje(16,n,m)
        ficje(17,n,m)=signy*ficje(17,n,m)
        ficje(18,n,m)=signz*ficje(18,n,m)
        ficje(19,n,m)=signx*ficje(19,n,m)
        ficje(20,n,m)=signy*ficje(20,n,m)
        ficje(21,n,m)=signz*ficje(21,n,m)
        ficje(25,n,m)=signx*ficje(25,n,m)
        ficje(28,n,m)=signy*ficje(28,n,m)
        ficje(30,n,m)=signz*ficje(30,n,m)
        ficje(31,n,m)=signx*ficje(31,n,m)
        ficje(34,n,m)=signy*ficje(34,n,m)
        ficje(36,n,m)=signz*ficje(36,n,m)
        ficje(37,n,m)=signx*ficje(37,n,m)
        ficje(40,n,m)=signy*ficje(40,n,m)
        ficje(42,n,m)=signz*ficje(42,n,m)
      enddo

      do n=0,nicjc(m)
        signx=sign(1.0,ficjc(10,n,m))
        signy=sign(1.0,ficjc(11,n,m))
        signz=sign(1.0,ficjc(12,n,m))
        ficjc(13,n,m)=signx*ficjc(13,n,m)
        ficjc(14,n,m)=signy*ficjc(14,n,m)
        ficjc(15,n,m)=signz*ficjc(15,n,m)
        ficjc(16,n,m)=signx*ficjc(16,n,m)
        ficjc(17,n,m)=signy*ficjc(17,n,m)
        ficjc(18,n,m)=signz*ficjc(18,n,m)
        ficjc(19,n,m)=signx*ficjc(19,n,m)
        ficjc(20,n,m)=signy*ficjc(20,n,m)
        ficjc(21,n,m)=signz*ficjc(21,n,m)
        ficjc(25,n,m)=signx*ficjc(25,n,m)
        ficjc(28,n,m)=signy*ficjc(28,n,m)
        ficjc(30,n,m)=signz*ficjc(30,n,m)
        ficjc(31,n,m)=signx*ficjc(31,n,m)
        ficjc(34,n,m)=signy*ficjc(34,n,m)
        ficjc(36,n,m)=signz*ficjc(36,n,m)
        ficjc(37,n,m)=signx*ficjc(37,n,m)
        ficjc(40,n,m)=signy*ficjc(40,n,m)
        ficjc(42,n,m)=signz*ficjc(42,n,m)
      enddo

      do n=0,niekc(m)
        signx=sign(1.0,fiekc(10,n,m))
        signy=sign(1.0,fiekc(11,n,m))
        signz=sign(1.0,fiekc(12,n,m))
        fiekc(13,n,m)=signx*fiekc(13,n,m)
        fiekc(14,n,m)=signy*fiekc(14,n,m)
        fiekc(15,n,m)=signz*fiekc(15,n,m)
        fiekc(16,n,m)=signx*fiekc(16,n,m)
        fiekc(17,n,m)=signy*fiekc(17,n,m)
        fiekc(18,n,m)=signz*fiekc(18,n,m)
        fiekc(19,n,m)=signx*fiekc(19,n,m)
        fiekc(20,n,m)=signy*fiekc(20,n,m)
        fiekc(21,n,m)=signz*fiekc(21,n,m)
        fiekc(25,n,m)=signx*fiekc(25,n,m)
        fiekc(28,n,m)=signy*fiekc(28,n,m)
        fiekc(30,n,m)=signz*fiekc(30,n,m)
        fiekc(31,n,m)=signx*fiekc(31,n,m)
        fiekc(34,n,m)=signy*fiekc(34,n,m)
        fiekc(36,n,m)=signz*fiekc(36,n,m)
        fiekc(37,n,m)=signx*fiekc(37,n,m)
        fiekc(40,n,m)=signy*fiekc(40,n,m)
        fiekc(42,n,m)=signz*fiekc(42,n,m)
      enddo

      do n=0,nicke(m)
        signx=sign(1.0,ficke(10,n,m))
        signy=sign(1.0,ficke(11,n,m))
        signz=sign(1.0,ficke(12,n,m))
        ficke(13,n,m)=signx*ficke(13,n,m)
        ficke(14,n,m)=signy*ficke(14,n,m)
        ficke(15,n,m)=signz*ficke(15,n,m)
        ficke(16,n,m)=signx*ficke(16,n,m)
        ficke(17,n,m)=signy*ficke(17,n,m)
        ficke(18,n,m)=signz*ficke(18,n,m)
        ficke(19,n,m)=signx*ficke(19,n,m)
        ficke(20,n,m)=signy*ficke(20,n,m)
        ficke(21,n,m)=signz*ficke(21,n,m)
        ficke(25,n,m)=signx*ficke(25,n,m)
        ficke(28,n,m)=signy*ficke(28,n,m)
        ficke(30,n,m)=signz*ficke(30,n,m)
        ficke(31,n,m)=signx*ficke(31,n,m)
        ficke(34,n,m)=signy*ficke(34,n,m)
        ficke(36,n,m)=signz*ficke(36,n,m)
        ficke(37,n,m)=signx*ficke(37,n,m)
        ficke(40,n,m)=signy*ficke(40,n,m)
        ficke(42,n,m)=signz*ficke(42,n,m)
      enddo

      do n=0,nickc(m)
        signx=sign(1.0,fickc(10,n,m))
        signy=sign(1.0,fickc(11,n,m))
        signz=sign(1.0,fickc(12,n,m))
        fickc(13,n,m)=signx*fickc(13,n,m)
        fickc(14,n,m)=signy*fickc(14,n,m)
        fickc(15,n,m)=signz*fickc(15,n,m)
        fickc(16,n,m)=signx*fickc(16,n,m)
        fickc(17,n,m)=signy*fickc(17,n,m)
        fickc(18,n,m)=signz*fickc(18,n,m)
        fickc(19,n,m)=signx*fickc(19,n,m)
        fickc(20,n,m)=signy*fickc(20,n,m)
        fickc(21,n,m)=signz*fickc(21,n,m)
        fickc(25,n,m)=signx*fickc(25,n,m)
        fickc(28,n,m)=signy*fickc(28,n,m)
        fickc(30,n,m)=signz*fickc(30,n,m)
        fickc(31,n,m)=signx*fickc(31,n,m)
        fickc(34,n,m)=signy*fickc(34,n,m)
        fickc(36,n,m)=signz*fickc(36,n,m)
        fickc(37,n,m)=signx*fickc(37,n,m)
        fickc(40,n,m)=signy*fickc(40,n,m)
        fickc(42,n,m)=signz*fickc(42,n,m)
      enddo

      do n=0,njekc(m)
        signx=sign(1.0,fjekc(10,n,m))
        signy=sign(1.0,fjekc(11,n,m))
        signz=sign(1.0,fjekc(12,n,m))
        fjekc(13,n,m)=signx*fjekc(13,n,m)
        fjekc(14,n,m)=signy*fjekc(14,n,m)
        fjekc(15,n,m)=signz*fjekc(15,n,m)
        fjekc(16,n,m)=signx*fjekc(16,n,m)
        fjekc(17,n,m)=signy*fjekc(17,n,m)
        fjekc(18,n,m)=signz*fjekc(18,n,m)
        fjekc(19,n,m)=signx*fjekc(19,n,m)
        fjekc(20,n,m)=signy*fjekc(20,n,m)
        fjekc(21,n,m)=signz*fjekc(21,n,m)
        fjekc(25,n,m)=signx*fjekc(25,n,m)
        fjekc(28,n,m)=signy*fjekc(28,n,m)
        fjekc(30,n,m)=signz*fjekc(30,n,m)
        fjekc(31,n,m)=signx*fjekc(31,n,m)
        fjekc(34,n,m)=signy*fjekc(34,n,m)
        fjekc(36,n,m)=signz*fjekc(36,n,m)
        fjekc(37,n,m)=signx*fjekc(37,n,m)
        fjekc(40,n,m)=signy*fjekc(40,n,m)
        fjekc(42,n,m)=signz*fjekc(42,n,m)
      enddo

      do n=0,njcke(m)
        signx=sign(1.0,fjcke(10,n,m))
        signy=sign(1.0,fjcke(11,n,m))
        signz=sign(1.0,fjcke(12,n,m))
        fjcke(14,n,m)=signy*fjcke(14,n,m)
        fjcke(15,n,m)=signz*fjcke(15,n,m)
        fjcke(16,n,m)=signx*fjcke(16,n,m)
        fjcke(17,n,m)=signy*fjcke(17,n,m)
        fjcke(18,n,m)=signz*fjcke(18,n,m)
        fjcke(19,n,m)=signx*fjcke(19,n,m)
        fjcke(20,n,m)=signy*fjcke(20,n,m)
        fjcke(21,n,m)=signz*fjcke(21,n,m)
        fjcke(25,n,m)=signx*fjcke(25,n,m)
        fjcke(28,n,m)=signy*fjcke(28,n,m)
        fjcke(30,n,m)=signz*fjcke(30,n,m)
        fjcke(31,n,m)=signx*fjcke(31,n,m)
        fjcke(34,n,m)=signy*fjcke(34,n,m)
        fjcke(36,n,m)=signz*fjcke(36,n,m)
        fjcke(37,n,m)=signx*fjcke(37,n,m)
        fjcke(40,n,m)=signy*fjcke(40,n,m)
        fjcke(42,n,m)=signz*fjcke(42,n,m)
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        signy=sign(1.0,fjckc(11,n,m))
        signz=sign(1.0,fjckc(12,n,m))
        fjckc(13,n,m)=signx*fjckc(13,n,m)
        fjckc(14,n,m)=signy*fjckc(14,n,m)
        fjckc(15,n,m)=signz*fjckc(15,n,m)
        fjckc(16,n,m)=signx*fjckc(16,n,m)
        fjckc(17,n,m)=signy*fjckc(17,n,m)
        fjckc(18,n,m)=signz*fjckc(18,n,m)
        fjckc(19,n,m)=signx*fjckc(19,n,m)
        fjckc(20,n,m)=signy*fjckc(20,n,m)
        fjckc(21,n,m)=signz*fjckc(21,n,m)
        fjckc(25,n,m)=signx*fjckc(25,n,m)
        fjckc(28,n,m)=signy*fjckc(28,n,m)
        fjckc(30,n,m)=signz*fjckc(30,n,m)
        fjckc(31,n,m)=signx*fjckc(31,n,m)
        fjckc(34,n,m)=signy*fjckc(34,n,m)
        fjckc(36,n,m)=signz*fjckc(36,n,m)
        fjckc(37,n,m)=signx*fjckc(37,n,m)
        fjckc(40,n,m)=signy*fjckc(40,n,m)
        fjckc(42,n,m)=signz*fjckc(42,n,m)
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------

         

