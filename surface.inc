c  surface.inc

      integer i01(ms)
      integer niejc(ms),nicje(ms),nicjc(ms)
      integer niekc(ms),nicke(ms),nickc(ms)
      integer njekc(ms),njcke(ms),njckc(ms)

      integer iejc(4,0:ij,ms),icje(4,0:ij,ms),icjc(4,0:ij,ms)
      integer iekc(4,0:ik,ms),icke(4,0:ik,ms),ickc(4,0:ik,ms)
      integer jekc(4,0:jk,ms),jcke(4,0:jk,ms),jckc(4,0:jk,ms)

      integer kuiejc(0:ij,ms),kvicje(0:ij,ms),kwicjc(0:ij,ms)
      integer kwiejc(0:ij,ms),kwicje(0:ij,ms),kuicje(0:ij,ms)
      integer juiekc(0:ik,ms),jwicke(0:ik,ms),jvickc(0:ik,ms)
      integer jviekc(0:ik,ms),jvicke(0:ik,ms),jwiekc(0:ik,ms)
      integer ivjekc(0:jk,ms),iwjcke(0:jk,ms),iujckc(0:jk,ms)
      integer iujekc(0:jk,ms),iujcke(0:jk,ms),ivjcke(0:jk,ms)

c     integer ijdex(nx,ny,nz),ikdex(nx,ny,nz),jkdex(nx,ny,nz)

      real*8 alfa10(0:ns1),alfa11(0:ns1),alfa20(0:ns2)
      real*8 sk0(ms),sk1(ms),sk2(ms),sk3(ms)
      real*8 ea(ms),eb(ms),ec(ms)
      real*8 phase(ms),phasep(ms),phasey(ms),phaser(ms)
      real*8 xsc(ms),ysc(ms),zsc(ms),pitch(ms),yaw(ms),roll(ms)
      real*8 xsc0(ms),ysc0(ms),zsc0(ms),pitch0(ms),yaw0(ms),roll0(ms)
      real*8 xsct(ms),ysct(ms),zsct(ms),pitcht(ms),yawt(ms),rollt(ms)
      real*8 xsctt(ms),ysctt(ms),zsctt(ms)

      real*8 xscn(ms),yscn(ms),zscn(ms),xsctn(ms),ysctn(ms),zsctn(ms)
      real*8 xscrk(4,ms),yscrk(4,ms),zscrk(4,ms)
      real*8 xsctrk(4,ms),ysctrk(4,ms),zsctrk(4,ms)

      real*8 xs(0:ns1,0:ns2,ms),ys(0:ns1,0:ns2,ms),zs(0:ns1,0:ns2,ms)
      real*8 xsn(0:ns1,0:ns2,ms),ysn(0:ns1,0:ns2,ms),zsn(0:ns1,0:ns2,ms)
      real*8 xse(0:ns1,0:ns2,ms),yse(0:ns1,0:ns2,ms),zse(0:ns1,0:ns2,ms)
      real*8 xss(0:ns1,0:ns2,ms),yss(0:ns1,0:ns2,ms),zss(0:ns1,0:ns2,ms)
      real*8 fst(0:ns1,0:ns2,12,ms),snd(0:ns1,0:ns2,15,ms)

      real*8 us(0:ns1,0:ns2,ms),vs(0:ns1,0:ns2,ms),ws(0:ns1,0:ns2,ms)
      real*8 use(0:ns1,0:ns2,ms),vse(0:ns1,0:ns2,ms),wse(0:ns1,0:ns2,ms)
      real*8 usrk(4,0:ns1,0:ns2,ms)
      real*8 vsrk(4,0:ns1,0:ns2,ms)
      real*8 wsrk(4,0:ns1,0:ns2,ms)

      real*8 fx(0:ns1,0:ns2,ms),fy(0:ns1,0:ns2,ms),fz(0:ns1,0:ns2,ms)
      real*8 fr(0:ns1,0:ns2,20,ms)
      real*8 fiejc(62,0:ij,ms),ficje(62,0:ij,ms),ficjc(76,0:ij,ms)
      real*8 fiekc(62,0:ik,ms),ficke(62,0:ik,ms),fickc(76,0:ik,ms)
      real*8 fjekc(62,0:jk,ms),fjcke(62,0:jk,ms),fjckc(76,0:jk,ms)

      real*8 ujc1(0:ns1,0:ns2,3,ms),ujc2(0:ns1,0:ns2,6,ms)
      real*8 vjc1(0:ns1,0:ns2,3,ms),vjc2(0:ns1,0:ns2,6,ms)
      real*8 wjc1(0:ns1,0:ns2,3,ms),wjc2(0:ns1,0:ns2,6,ms)
      real*8 pjc0(0:ns1,0:ns2,ms),pjc(0:ns1,0:ns2,6,ms)
      real*8 pjc1(0:ns1,0:ns2,3,ms),pjc2(0:ns1,0:ns2,6,ms)

      real*8 ukiejc(0:ij,ms),vkicje(0:ij,ms),wkicjc(0:ij,ms)
      real*8 wkiejc(0:ij,ms),wkicje(0:ij,ms),ukicje(0:ij,ms)
      real*8 ujiekc(0:ik,ms),wjicke(0:ik,ms),vjickc(0:ik,ms)
      real*8 vjiekc(0:ik,ms),vjicke(0:ik,ms),wjiekc(0:ik,ms)
      real*8 vijekc(0:jk,ms),wijcke(0:jk,ms),uijckc(0:jk,ms)
      real*8 uijekc(0:jk,ms),uijcke(0:jk,ms),vijcke(0:jk,ms)

      real*8 uudx(0:jk,ms),uvdy(0:ik,ms),uwdz(0:ij,ms),pdx(0:jk,ms)
      real*8 udxx(2,0:jk,ms),udyy(2,0:ik,ms),udzz(2,0:ij,ms)
      real*8 vudx(0:jk,ms),vvdy(0:ik,ms),vwdz(0:ij,ms),pdy(0:ik,ms)
      real*8 vdxx(2,0:jk,ms),vdyy(2,0:ik,ms),vdzz(2,0:ij,ms)
      real*8 wudx(0:jk,ms),wvdy(0:ik,ms),wwdz(0:ij,ms),pdz(0:ij,ms)
      real*8 wdxx(2,0:jk,ms),wdyy(2,0:ik,ms),wdzz(2,0:ij,ms)
      real*8 udx(0:jk,ms),vdx(0:jk,ms),wdx(0:jk,ms)
      real*8 udy(0:ik,ms),vdy(0:ik,ms),wdy(0:ik,ms)
      real*8 udz(0:ij,ms),vdz(0:ij,ms),wdz(0:ij,ms)
      real*8 pdxx(2,0:jk,ms),pdyy(2,0:ik,ms),pdzz(2,0:ij,ms)

      common/ncount/i01,niejc,nicje,nicjc,
     .              niekc,nicke,nickc,njekc,njcke,njckc
      common/index/iejc,icje,icjc,iekc,icke,ickc,jekc,jcke,jckc,
     .             kuiejc,kvicje,kwicjc,kwiejc,kwicje,kuicje,
     .             juiekc,jwicke,jvickc,jviekc,jvicke,jwiekc,
     .             ivjekc,iwjcke,iujckc,iujekc,iujcke,ivjcke
c    .             ,ijdex,ikdex,jkdex
      common/surface/ea,eb,ec,phase,phasep,phasey,phaser,
     .               xsc,ysc,zsc,pitch,yaw,roll,
     .               xsc0,ysc0,zsc0,pitch0,yaw0,roll0,
     .               xsct,ysct,zsct,pitcht,yawt,rollt,xsctt,ysctt,zsctt,
     .               alfa10,alfa11,alfa20,xs,ys,zs,xsn,ysn,zsn,
     .               xse,yse,zse,xss,yss,zss,fst,snd
      common/fijk/fiejc,ficje,ficjc,fiekc,ficke,fickc,fjekc,fjcke,fjckc
      common/sf/sk0,sk1,sk2,sk3,fx,fy,fz,fr
      common/suvw/us,vs,ws,use,vse,wse,usrk,vsrk,wsrk
      common/prescribe/xscn,yscn,zscn,xsctn,ysctn,zsctn,
     .                 xscrk,yscrk,zscrk,xsctrk,ysctrk,zsctrk
      common/jc/ujc1,vjc1,wjc1,ujc2,vjc2,wjc2,pjc0,pjc1,pjc2,pjc
      common/ci/ukiejc,vkicje,wkicjc,wkiejc,wkicje,ukicje,
     .          ujiekc,wjicke,vjickc,vjiekc,vjicke,wjiekc,
     .          vijekc,wijcke,uijckc,uijekc,uijcke,vijcke
      common/cd/uudx,uvdy,uwdz,pdx,udxx,udyy,udzz,
     .          vudx,vvdy,vwdz,pdy,vdxx,vdyy,vdzz,
     .          wudx,wvdy,wwdz,pdz,wdxx,wdyy,wdzz,
     .          udx,vdx,wdx,udy,vdy,wdy,udz,vdz,wdz,pdxx,pdyy,pdzz
