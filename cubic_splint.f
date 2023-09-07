c-----------------------------------------------------------------------
c
      subroutine cubic_splint(xx,yy,cc,ns,nn,num)
      integer ns,nn,num
      real*8 dm,dn
      real*8 xx(0:nn),yy(0:nn,num),cc(0:nn,num)
      real*8 h(0:ns-1),dyy(0:ns-1,num)
      real*8 a(ns),b(ns),d(ns,num)

c  coefficient matrix
c
      do m=0,ns-1
        h(m)=xx(m+1)-xx(m)
      enddo

      do m=1,ns-1
        a(m)=h(m-1)/(h(m)+h(m-1))
        b(m)=1.0-a(m)
      enddo
      a(ns)=h(0)/(h(0)+h(ns-1))
      b(ns)=1.0d0-a(ns)

      do n=1,num

        do m=0,ns-1     
          dyy(m,n)=yy(m+1,n)-yy(m,n)
        enddo

        do m=1,ns-1
          d(m,n)=6.0d0*(dyy(m,n)/h(m)-dyy(m-1,n)/h(m-1))/(h(m)+h(m-1))
        enddo
        d(ns,n)=6.0d0*(dyy(0,n)/h(0)-dyy(ns-1,n)/h(ns-1))/(h(0)+h(ns-1))

      enddo

      a(1)=a(1)/2.0d0
      b(1)=b(1)/2.0d0
      do n=1,num
        d(1,n)=d(1,n)/2.0d0
      enddo
      b(ns)=b(ns)/2.0d0
      a(ns)=a(ns)/2.0d0
      do n=1,num
        d(ns,n)=d(ns,n)/2.0d0
      enddo
 
c  elimination to upper triangular matrix
c
c      
c |  1                                        -1       |
c | a(1)    1     b(1)                                 |
c |  0     a(2)    2      b(2)      		       |
c |						       |
c |                      ...		    	       |
c |						       |
c |                             a(ns-1)   2   b(ns-1)  |
c |  0     b(ns)                         a(ns)   1     |
c
c
c  ======>
c
c
c |  1                                        -1       |
c |  0      1     b(1)                       a(1)      |
c |  0      0       1      b(2)              a(2)      |
c |                                                    |
c |                      ...                           |
c |                               1   b(ns-2) a(ns-2)  |
c |                            a(ns-1)   2    b(ns-1)  |
c |  0                          b(ns)   a(ns)    1     |
c
c
c  ======>
c
c
c |  1                                        -1       |
c |  0      1     b(1)                       a(1)      |
c |  0      0       1      b(2)              a(2)      |
c |                                                    |
c |                      ...                           |
c |                               1   b(ns-2) a(ns-2)  |
c |                               0     1    a(ns-1)   |
c |  0                            0     0      1       |
c
c
      do m=2,ns-2
        dm=2.0d0-b(m-1)*a(m)
        do n=1,num
          d(m,n)=(d(m,n)-d(m-1,n)*a(m))/dm
        enddo
        a(m)=-a(m-1)*a(m)/dm
        b(m)=b(m)/dm

        dn=1.0d0-a(m-1)*b(ns)
        do n=1,num
          d(ns,n)=(d(ns,n)-d(m-1,n)*b(ns))/dn
        enddo
        b(ns)=-b(m-1)*b(ns)/dn
        a(ns)=a(ns)/dn
      enddo

      dm=2.0d0-b(ns-2)*a(ns-1)
      do n=1,num
        d(ns-1,n)=(d(ns-1,n)-d(ns-2,n)*a(ns-1))/dm
      enddo
      a(ns-1)=(b(ns-1)-a(ns-2)*a(ns-1))/dm

      dn=1.0d0-a(ns-2)*b(ns)
      do n=1,num
        d(ns,n)=(d(ns,n)-d(ns-2,n)*b(ns))/dn
      enddo
      a(ns)=(a(ns)-b(ns-2)*b(ns))/dn

      dn=1.0d0-a(ns-1)*a(ns)
      do n=1,num
        d(ns,n)=(d(ns,n)-d(ns-1,n)*a(ns))/dn
      enddo

c  back substitution
c
      do n=1,num
        cc(ns,n)=d(ns,n)
        cc(ns-1,n)=d(ns-1,n)-a(ns-1)*cc(ns,n)
        do m=ns-2,1,-1
          cc(m,n)=d(m,n)-a(m)*cc(ns,n)-b(m)*cc(m+1,n)
        enddo
        cc(0,n)=cc(ns,n)
      enddo
 
      return
      end


c-----------------------------------------------------------------------
