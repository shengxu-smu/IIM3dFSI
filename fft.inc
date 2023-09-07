c  fft.inc
      real*8 wsave1(ns1+15),wsave2(ns2+15),wsave3(2*ns1+2+15)
      real*8 r1(ns2,ns1),w1(ns2,ns1)
      real*8 r2(ns1,ns2),w2(ns1,ns2)
      real*8 r3(ns2/2,2*ns1+2),w3(ns2/2,2*ns1+2)
      real*8 r4(ns1+1,ns2),w4(ns1+1,ns2)

      real*8 wsavei(3*nx+15),wsavej(3*ny+15),wsavek(3*nz+15)
      real*8 ri(ny,nx),wi(ny,nx)
      real*8 rj(nx,ny),wj(nx,ny)
      real*8 rk(ny,nz),wk(ny,nz)

      common/fft/wsave1,wsave2,wsave3,r1,w1,r2,w2,r3,w3,r4,w4,
     .           wsavei,wsavej,wsavek,ri,wi,rj,wj,rk,wk

