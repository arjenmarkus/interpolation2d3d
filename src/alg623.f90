module alg623
    implicit none
contains
!*==ADNODE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
!     ALGORITHM 623 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 4,
!     DEC., 1984, P. 437.
      subroutine adnode(kk,x,y,z,iadj,iend,ier)
!*--********************************************************************
!A INPUT  - KK
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - IADJ
!A INPUT  - IEND
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       BDYADD   COVSPH   INDEX    INTADD   SWAP     SWPTST   TRFIND
! called by   TRMESH   TRMESH
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DUM      I1       I2       I3       IN1      IND21    IND2F    INDK1    INDKF    INDKL    IO1      IO2      K        KM1      NABOR1   P
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      !integer index
!*** End of declarations inserted by SPAG
      integer kk , iadj(*) , iend(kk) , ier
      real x(kk) , y(kk) , z(kk)
      !logical swptst
      !external index
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS NODE KK TO A TRIANGULATION OF THE
! CONVEX HULL OF NODES 1,...,KK-1, PRODUCING A TRIANGULATION
! OF THE CONVEX HULL OF NODES 1,...,KK.  A SEQUENCE OF EDGE
! SWAPS IS THEN APPLIED TO THE MESH, RESULTING IN AN OPTIMAL
! TRIANGULATION.  ADNODE IS PART OF AN INTERPOLATION PACKAGE
! WHICH ALSO PROVIDES ROUTINES TO INITIALIZE THE DATA STRUC-
! TURE, PLOT THE MESH, AND DELETE ARCS.
!
! INPUT PARAMETERS -   KK - INDEX OF THE NODE TO BE ADDED
!                           TO THE MESH.  KK .GE. 4.
!
!                   X,Y,Z - VECTORS OF LENGTH .GE. KK CON-
!                           TAINING CARTESIAN COORDINATES
!                           OF THE NODES.  (X(I),Y(I),Z(I))
!                           DEFINES NODE I FOR I = 1,...,KK.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           1,...,KK-1.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS IN IADJ FOR
!                           EACH NODE IN THE MESH.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! KK, X, Y, AND Z ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.
!
!                           IER - ERROR INDICATOR
!                                 IER = 0 IF NO ERRORS
!                                         WERE ENCOUNTERED.
!                                 IER = 1 IF ALL NODES
!                                         (INCLUDING KK) ARE
!                                         COLLINEAR.
!
! MODULES REFERENCED BY ADNODE - TRFIND, INTADD, BDYADD,
!                                COVSPH, SHIFTD, INDEX,
!                                SWPTST, SWAP
!
!***********************************************************
!
      integer k , km1 , i1 , i2 , i3 , indkf , indkl , nabor1 , io1 ,   &
     &        io2 , in1 , indk1 , ind2f , ind21
      real p(3) , dum
!
! LOCAL PARAMETERS -
!
! K =        LOCAL COPY OF KK
! KM1 =      K - 1
! I1,I2,I3 = VERTICES OF A TRIANGLE CONTAINING K
! INDKF =    IADJ INDEX OF THE FIRST NEIGHBOR OF K
! INDKL =    IADJ INDEX OF THE LAST NEIGHBOR OF K
! NABOR1 =   FIRST NEIGHBOR OF K BEFORE ANY SWAPS OCCUR
! IO1,IO2 =  ADJACENT NEIGHBORS OF K DEFINING AN ARC TO
!              BE TESTED FOR A SWAP
! IN1 =      VERTEX OPPOSITE K -- FIRST NEIGHBOR OF IO2
!              WHICH PRECEDES IO1.  IN1,IO1,IO2 ARE IN
!              COUNTERCLOCKWISE ORDER.
! INDK1 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR K
! IND2F =    INDEX OF THE FIRST NEIGHBOR OF IO2
! IND21 =    INDEX OF IO1 IN THE ADJACENCY LIST FOR IO2
! P =        CARTESIAN COORDINATES OF NODE KK
! DUM =      DUMMY PARAMETER FOR CALL TO TRFIND
!
      ier = 0
      k = kk
!
! INITIALIZATION
!
      km1 = k - 1
      p(1) = x(k)
      p(2) = y(k)
      p(3) = z(k)
!
! ADD NODE K TO THE MESH
!
      call trfind(km1,p,x,y,z,iadj,iend,dum,dum,dum,i1,i2,i3)
      if ( i1==0 ) then
!
! ALL NODES ARE COLLINEAR
!
         ier = 1
         goto 99999
      else
         if ( i3/=0 ) call intadd(k,i1,i2,i3,iadj,iend)
         if ( i3==0 .and. i1/=i2 ) call bdyadd(k,i1,i2,iadj,iend)
         if ( i1==i2 ) call covsph(k,i1,iadj,iend)
!
! INITIALIZE VARIABLES FOR OPTIMIZATION OF THE MESH
!
         indkf = iend(km1) + 1
         indkl = iend(k)
         nabor1 = iadj(indkf)
         io2 = nabor1
         indk1 = indkf + 1
         io1 = iadj(indk1)
      endif
!
! BEGIN LOOP -- FIND THE VERTEX OPPOSITE K
!
 100  ind2f = 1
      if ( io2/=1 ) ind2f = iend(io2-1) + 1
      ind21 = indexr(io2,io1,iadj,iend)
      if ( ind2f==ind21 ) then
!
! IN1 IS THE LAST NEIGHBOR OF IO2
!
         ind21 = iend(io2)
         in1 = iadj(ind21)
         if ( in1==0 ) goto 200
      else
         in1 = iadj(ind21-1)
      endif
!
! SWAP TEST -- IF A SWAP OCCURS, TWO NEW ARCS ARE OPPOSITE K
!              AND MUST BE TESTED.  INDK1 AND INDKF MUST BE
!              DECREMENTED.
!
      if ( swptst(io1,io2,in1,k,x,y,z) ) then
         call swap(in1,k,io1,io2,iadj,iend)
         io1 = in1
         indk1 = indk1 - 1
         indkf = indkf - 1
         goto 100
      endif
!
! NO SWAP OCCURRED.  RESET IO2 AND IO1, AND TEST FOR
!   TERMINATION.
!
 200  if ( io1==nabor1 ) return
      io2 = io1
      indk1 = indk1 + 1
      if ( indk1>indkl ) indk1 = indkf
      io1 = iadj(indk1)
      if ( io1/=0 ) goto 100
      return
99999 end subroutine adnode
!*==APLYR.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine aplyr(x,y,z,cx,sx,cy,sy,xp,yp,zp)
!*--********************************************************************
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - CX
!A INPUT  - SX
!A INPUT  - CY
!A INPUT  - SY
!A OUTPUT - XP
!A OUTPUT - YP
!A OUTPUT - ZP
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  T
! uses PARAMs *** NONE ****
!*++********************************************************************
      real x , y , z , cx , sx , cy , sy , xp , yp , zp
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE ROTATION R DEFINED BY SUB-
! ROUTINE CONSTR TO THE UNIT VECTOR (X Y Z)**T, I.E. (X,Y,Z)
! IS ROTATED TO (XP,YP,ZP).  IF (XP,YP,ZP) LIES IN THE
! SOUTHERN HEMISPHERE (ZP .LT. 0), (XP,YP) ARE SET TO THE
! COORDINATES OF THE NEAREST POINT OF THE EQUATOR, ZP RE-
! MAINING UNCHANGED.
!
! INPUT PARAMETERS - X,Y,Z - COORDINATES OF A POINT ON THE
!                            UNIT SPHERE.
!
!              CX,SX,CY,SY - ELEMENTS OF THE ROTATION DE-
!                            FINED BY CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED EXCEPT AS NOTED BELOW.
!
! OUTPUT PARAMETERS - XP,YP,ZP - COORDINATES OF THE ROTATED
!                                POINT ON THE SPHERE UNLESS
!                                ZP .LT. 0, IN WHICH CASE
!                                (XP,YP,0) IS THE CLOSEST
!                                POINT OF THE EQUATOR TO THE
!                                ROTATED POINT.  STORAGE FOR
!                                XP, YP, AND ZP MAY COINCIDE
!                                WITH STORAGE FOR X, Y, AND
!                                Z, RESPECTIVELY, IF THE
!                                LATTER NEED NOT BE SAVED.
!
! MODULES REFERENCED BY APLYR - NONE
!
! INTRINSIC FUNCTION CALLED BY APLYR - SQRT
!
!***********************************************************
!
      real t
!
! LOCAL PARAMETER -
!
! T = TEMPORARY VARIABLE
!
      t = sx*y + cx*z
      yp = cx*y - sx*z
      zp = sy*x + cy*t
      xp = cy*x - sy*t
      if ( zp>=0. ) return
!
! MOVE (XP,YP,ZP) TO THE EQUATOR
!
      t = sqrt(xp*xp+yp*yp)
      if ( t==0. ) then
!
! MOVE THE SOUTH POLE TO AN ARBITRARY POINT OF THE EQUATOR
!
         xp = 1.
         yp = 0.
         goto 99999
      endif
      xp = xp/t
      yp = yp/t
      return
99999 end subroutine aplyr
!*==APLYRT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine aplyrt(g1p,g2p,cx,sx,cy,sy,g)
!*--********************************************************************
!A INPUT  - G1P
!A INPUT  - G2P
!A INPUT  - CX
!A INPUT  - SX
!A INPUT  - CY
!A INPUT  - SY
!A OUTPUT - G
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADG    GRADL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  T
! uses PARAMs *** NONE ****
!*++********************************************************************
      real g1p , g2p , cx , sx , cy , sy , g(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE APPLIES THE INVERSE (TRANSPOSE) OF THE
! ROTATION DEFINED BY SUBROUTINE CONSTR TO THE VECTOR
! (G1P G2P 0)**T, I.E. THE GRADIENT (G1P,G2P,0) IN THE ROT-
! ATED COORDINATE SYSTEM IS MAPPED TO (G1,G2,G3) IN THE
! ORIGINAL COORDINATE SYSTEM.
!
! INPUT PARAMETERS - G1P,G2P - X- AND Y-COMPONENTS, RESPECT-
!                              IVELY, OF THE GRADIENT IN THE
!                              ROTATED COORDINATE SYSTEM.
!
!                CX,SX,CY,SY - ELEMENTS OF THE ROTATION R
!                              CONSTRUCTED BY SUBROUTINE
!                              CONSTR.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - G - X-, Y-, AND Z-COMPONENTS (IN THAT
!                         ORDER) OF THE INVERSE ROTATION
!                         APPLIED TO (G1P,G2P,0) -- GRADIENT
!                         IN THE ORIGINAL COORDINATE SYSTEM.
!
! MODULES REFERENCED BY APLYRT - NONE
!
!***********************************************************
!
      real t
!
! LOCAL PARAMETERS -
!
! T = TEMPORARY VARIABLE
!
      t = sy*g1p
      g(1) = cy*g1p
      g(2) = cx*g2p - sx*t
      g(3) = -sx*g2p - cx*t
      end subroutine aplyrt
!*==ARCINT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine arcint(p,p1,p2,w1,w2,g1,g2,w,g,gn)
!*--********************************************************************
!A INPUT  - P
!A INPUT  - P1
!A INPUT  - P2
!A INPUT  - W1
!A INPUT  - W2
!A INPUT  - G1
!A INPUT  - G2
!A OUTPUT - W
!A OUTPUT - G
!A OUTPUT - GN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ARCLEN
! called by   INTRC1   WVAL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  A        AL       GT       I        LUN      S        T        TAU1     TAU2     UN       UNORM
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      !real arclen
!*** End of declarations inserted by SPAG
      real p(3) , p1(3) , p2(3) , w1 , w2 , g1(3) , g2(3) , w , g(3) ,  &
     &     gn
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN 3 POINTS P, P1, AND P2 LYING ON A COMMON GEODESIC
! OF THE UNIT SPHERE WITH P BETWEEN P1 AND P2, ALONG WITH
! DATA VALUES AND GRADIENTS AT P1 AND P2, THIS SUBROUTINE
! COMPUTES AN INTERPOLATED VALUE W AND A GRADIENT VECTOR G
! AT P.  W IS COMPUTED BY HERMITE CUBIC INTERPOLATION REL-
! ATIVE TO ARC-LENGTH ALONG THE GEODESIC.  THE TANGENTIAL
! COMPONENT OF G IS THE DERIVATIVE (WITH RESPECT TO ARC-
! LENGTH) OF THE CUBIC INTERPOLANT AT P, WHILE THE NORMAL
! COMPONENT OF G IS OBTAINED BY LINEAR INTERPOLATION OF THE
! NORMAL COMPONENTS OF THE GRADIENTS AT P1 AND P2.  THIS
! ALGORITHM IS DUE TO C. L. LAWSON.
!
! INPUT PARAMETERS - P - CARTESIAN COORDINATES OF A POINT
!                        LYING ON THE ARC DEFINED BY P1 AND
!                        P2.
!
!                P1,P2 - COORDINATES OF DISTINCT POINTS ON
!                        THE UNIT SPHERE DEFINING AN ARC
!                        WITH LENGTH LESS THAN 180 DEGREES.
!
!                W1,W2 - DATA VALUES ASSOCIATED WITH P1 AND
!                        P2, RESPECTIVELY.
!
!                G1,G2 - GRADIENT VECTORS ASSOCIATED WITH P1
!                        AND P2.  G1 AND G2 ARE ORTHOGONAL
!                        TO P1 AND P2, RESPECTIVELY.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                    G - ARRAY OF LENGTH 3.
!
! OUTPUT PARAMETERS - W - INTERPOLATED VALUE AT P.
!
!                     G - INTERPOLATED GRADIENT AT P.
!
!                    GN - NORMAL COMPONENT OF G WITH THE
!                         DIRECTION P1 X P2 TAKEN TO BE
!                         POSITIVE.  THE EXTRAPOLATION
!                         PROCEDURE REQUIRES THIS COMPONENT.
!
! FOR EACH VECTOR V, V(1), V(2), AND V(3) CONTAIN X-, Y-,
!   AND Z-COMPONENTS, RESPECTIVELY.
!
! MODULES REFERENCED BY ARCINT - ARCLEN
!
! INTRINSIC FUNCTION CALLED BY ARCINT - SQRT
!
!***********************************************************
!
      integer i , lun
      real un(3) , unorm , tau1 , tau2 , a , al , s , t , gt
      data lun/6/
!
! LOCAL PARAMETERS -
!
! I =         DO-LOOP INDEX
! LUN =       LOGICAL UNIT FOR ERROR MESSAGES
! UN =        UNIT NORMAL TO THE PLANE OF P, P1, AND P2
! UNORM =     EUCLIDEAN NORM OF P1 X P2 -- USED TO NORMALIZE
!               UN
! TAU1,TAU2 = TANGENTIAL DERIVATIVES (COMPONENTS OF G1,G2)
!               AT P1 AND P2
! A =         ANGLE IN RADIANS (ARC-LENGTH) BETWEEN P1 AND
!               P2
! AL =        ARC-LENGTH BETWEEN P1 AND P
! S =         NORMALIZED VALUE OF AL -- AS P VARIES FROM P1
!               TO P2, S VARIES FROM 0 TO 1
! T =         1-S -- S AND T ARE BARYCENTRIC COORDINATES OF
!               P WITH RESPECT TO THE ARC FROM P1 TO P2
! GT =        TANGENTIAL COMPONENT OF G -- COMPONENT IN THE
!               DIRECTION UN X P
!
!
! COMPUTE UNIT NORMAL UN
!
      un(1) = p1(2)*p2(3) - p1(3)*p2(2)
      un(2) = p1(3)*p2(1) - p1(1)*p2(3)
      un(3) = p1(1)*p2(2) - p1(2)*p2(1)
      unorm = sqrt(un(1)*un(1)+un(2)*un(2)+un(3)*un(3))
      if ( unorm/=0. ) then
!
! NORMALIZE UN
!
         do i = 1 , 3
            un(i) = un(i)/unorm
         enddo
!
! COMPUTE TANGENTIAL DERIVATIVES AT THE ENDPOINTS --
!   TAU1 = (G1,UN X P1) = (G1,P2)/UNORM AND
!   TAU2 = (G2,UN X P2) = -(G2,P1)/UNORM.
!
         tau1 = (g1(1)*p2(1)+g1(2)*p2(2)+g1(3)*p2(3))/unorm
         tau2 = -(g2(1)*p1(1)+g2(2)*p1(2)+g2(3)*p1(3))/unorm
!
! COMPUTE ARC-LENGTHS A, AL
!
         a = arclen(p1,p2)
         if ( a/=0. ) then
            al = arclen(p1,p)
!
! COMPUTE W BY HERMITE CUBIC INTERPOLATION
!
            s = al/a
            t = 1. - s
            w = w1*(2.*s+1.)*t*t + w2*(3.-2.*s)                         &
     &          *s*s + a*s*t*(tau1*t-tau2*s)
!
! COMPUTE TANGENTIAL AND NORMAL DERIVATIVES AT P
!
            gt = 6.*s*t*(w2-w1)/a + tau1*t*(1.-3.*s) + tau2*s*(3.*s-2.)
            gn = t*(un(1)*g1(1)+un(2)*g1(2)+un(3)*g1(3))                &
     &           + s*(un(1)*g2(1)+un(2)*g2(2)+un(3)*g2(3))
!
! COMPUTE G = GT*(UN X P) + GN*UN
!
            g(1) = gt*(un(2)*p(3)-un(3)*p(2)) + gn*un(1)
            g(2) = gt*(un(3)*p(1)-un(1)*p(3)) + gn*un(2)
            g(3) = gt*(un(1)*p(2)-un(2)*p(1)) + gn*un(3)
            return
         endif
      endif
!
! P1 X P2 = 0.  PRINT AN ERROR MESSAGE AND TERMINATE
!   PROCESSING.
!
      write (lun,99001) (p1(i),i=1,3) , (p2(i),i=1,3)
99001 format ('1','ERROR IN ARCINT -- P1 = ',2(f9.6,',  '),f9.6/' ',19x,&
     &        'P2 = ',2(f9.6,',  '),f9.6)
      stop
      end subroutine arcint
!*==ARCLEN.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      real function arclen(p,q)
!*--********************************************************************
!A INPUT  - P
!A INPUT  - Q
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   ARCINT   INTRC1
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  D        I
! uses PARAMs *** NONE ****
!*++********************************************************************
      real p(3) , q(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION COMPUTES THE ARC-LENGTH (ANGLE IN RADIANS)
! BETWEEN A PAIR OF POINTS ON THE UNIT SPHERE.
!
! INPUT PARAMETERS - P,Q - VECTORS OF LENGTH 3 CONTAINING
!                          THE X-, Y-, AND Z-COORDINATES (IN
!                          THAT ORDER) OF POINTS ON THE UNIT
!                          SPHERE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER - ARCLEN - ANGLE IN RADIANS BETWEEN THE
!                             UNIT VECTORS P AND Q.  0 .LE.
!                             ARCLEN .LE. PI.
!
! MODULES REFERENCED BY ARCLEN - NONE
!
! INTRINSIC FUNCTIONS CALLED BY ARCLEN - ATAN, SQRT
!
!***********************************************************
!
      integer i
      real d
!
! LOCAL PARAMETERS -
!
! I = DO-LOOP INDEX
! D = EUCLIDEAN NORM SQUARED OF P+Q
!
      d = 0.
      do i = 1 , 3
         d = d + (p(i)+q(i))**2
      enddo
      if ( d==0. ) then
!
! P AND Q ARE SEPARATED BY 180 DEGREES
!
         arclen = 4.*atan(1.)
         return
      elseif ( d>=4. ) then
!
! P AND Q COINCIDE
!
         arclen = 0.
         goto 99999
      endif
      arclen = 2.*atan(sqrt((4.-d)/d))
      return
99999 end function arclen
!*==BDYADD.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine bdyadd(kk,i1,i2,iadj,iend)
!*--********************************************************************
!A INPUT  - KK
!A INPUT  - I1
!A INPUT  - I2
!A OUTPUT - IADJ
!A OUTPUT - IEND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       SHIFTD
! called by   ADNODE   ADNODE
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IMAX     IMIN     INDX     K        KEND     KM1      N1       N2       NEXT     NF       NL       NLEFT    NRIGHT
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer kk , i1 , i2 , iadj(*) , iend(kk)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS A BOUNDARY NODE TO A TRIANGULATION OF
! A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
! ARE UPDATED WITH THE INSERTION OF NODE KK.
!
! INPUT PARAMETERS -   KK - INDEX OF AN EXTERIOR NODE TO BE
!                           ADDED.  KK .GE. 4.
!
!                      I1 - FIRST (RIGHTMOST AS VIEWED FROM
!                           KK) BOUNDARY NODE IN THE MESH
!                           WHICH IS VISIBLE FROM KK - THE
!                           LINE SEGMENT KK-I1 INTERSECTS
!                           NO ARCS.
!
!                      I2 - LAST (LEFTMOST) BOUNDARY NODE
!                           WHICH IS VISIBLE FROM KK.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           IN THE MESH.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS IN IADJ FOR
!                           EACH NODE IN THE MESH.
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
! THE VERTICES I1 AND I2.  I1 AND I2 MAY BE DETERMINED BY
! TRFIND.
!
! KK, I1, AND I2 ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  NODE KK WILL BE
!                                 CONNECTED TO I1, I2, AND
!                                 ALL BOUNDARY NODES BETWEEN
!                                 THEM.  NO OPTIMIZATION OF
!                                 THE MESH IS PERFORMED.
!
! MODULE REFERENCED BY BDYADD - SHIFTD
!
! INTRINSIC FUNCTIONS CALLED BY BDYADD - MIN0, MAX0
!
!***********************************************************
!
      integer k , km1 , nright , nleft , nf , nl , n1 , n2 , i , imin , &
     &        imax , kend , next , indx
!
! LOCAL PARAMETERS -
!
! K =            LOCAL COPY OF KK
! KM1 =          K - 1
! NRIGHT,NLEFT = LOCAL COPIES OF I1, I2
! NF,NL =        INDICES OF IADJ BOUNDING THE PORTION OF THE
!                  ARRAY TO BE SHIFTED
! N1 =           IADJ INDEX OF THE FIRST NEIGHBOR OF NLEFT
! N2 =           IADJ INDEX OF THE LAST NEIGHBOR OF NRIGHT
! I =            DO-LOOP INDEX
! IMIN,IMAX =    BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
!                  ELEMENTS OF IEND TO BE INCREMENTED
! KEND =         POINTER TO THE LAST NEIGHBOR OF K IN IADJ
! NEXT =         NEXT BOUNDARY NODE TO BE CONNECTED TO KK
! INDX =         INDEX FOR IADJ
!
      k = kk
      km1 = k - 1
      nright = i1
      nleft = i2
!
! INITIALIZE VARIABLES
!
      nl = iend(km1)
      n1 = 1
      if ( nleft/=1 ) n1 = iend(nleft-1) + 1
      n2 = iend(nright)
      nf = max0(n1,n2)
!
! INSERT K AS A NEIGHBOR OF MAX(NRIGHT,NLEFT)
!
      call shiftd(nf,nl,2,iadj)
      iadj(nf+1) = k
      imin = max0(nright,nleft)
      do i = imin , km1
         iend(i) = iend(i) + 2
      enddo
!
! INITIALIZE KEND AND INSERT K AS A NEIGHBOR OF
!   MIN(NRIGHT,NLEFT)
!
      kend = nl + 3
      nl = nf - 1
      nf = min0(n1,n2)
      call shiftd(nf,nl,1,iadj)
      iadj(nf) = k
      imax = imin - 1
      imin = min0(nright,nleft)
      do i = imin , imax
         iend(i) = iend(i) + 1
      enddo
!
! INSERT NRIGHT AS THE FIRST NEIGHBOR OF K
!
      iadj(kend) = nright
!
! INITIALIZE INDX FOR LOOP ON BOUNDARY NODES BETWEEN NRIGHT
!   AND NLEFT
!
      indx = iend(nright) - 2
      do
         next = iadj(indx)
         if ( next==nleft ) then
!
! INSERT NLEFT AND 0 AS THE LAST NEIGHBORS OF K
!
            iadj(kend+1) = nleft
            kend = kend + 2
            iadj(kend) = 0
            iend(k) = kend
            exit
         else
!
! CONNECT NEXT AND K
!
            kend = kend + 1
            iadj(kend) = next
            indx = iend(next)
            iadj(indx) = k
            indx = indx - 1
         endif
      enddo
      end subroutine bdyadd
!*==BNODES.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine bnodes(n,iadj,iend,nb,na,nt,nodes)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - IADJ
!A INPUT  - IEND
!A OUTPUT - NB
!A OUTPUT - NA
!A OUTPUT - NT
!A OUTPUT - NODES
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0001   AA0003
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  INDF     INDL     K        N0       NN       NST
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iadj(*) , iend(n) , nb , na , nt , nodes(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N POINTS ON THE UNIT SPHERE,
! THIS SUBROUTINE RETURNS A VECTOR CONTAINING THE INDICES
! (IF ANY) OF THE COUNTERCLOCKWISE-ORDERED SEQUENCE OF NODES
! ON THE BOUNDARY OF THE CONVEX HULL OF THE SET OF POINTS.
! THE BOUNDARY IS EMPTY IF THE POINTS DO NOT LIE IN A SINGLE
! HEMISPHERE.  THE NUMBERS OF BOUNDARY NODES, ARCS, AND
! TRIANGLES ARE ALSO RETURNED.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!
!                     IADJ - SET OF ADJACENCY LISTS OF
!                            NODES IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                    NODES - VECTOR OF LENGTH .GE. NB.
!                            (NB .LE. N).
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND ARE NOT
! ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -   NB - NUMBER OF BOUNDARY NODES.
!
!                    NA,NT - NUMBER OF ARCS AND TRIANGLES,
!                            RESPECTIVELY, IN THE MESH.
!
!                    NODES - VECTOR OF NB BOUNDARY NODE
!                            INDICES RANGING FROM 1 TO N.
!
! MODULES REFERENCED BY BNODES - NONE
!
!***********************************************************
!
      integer nn , nst , indl , k , n0 , indf
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! NST =  FIRST ELEMENT OF NODES -- ARBITRARILY CHOSEN
! INDL = IADJ INDEX OF THE LAST NEIGHBOR OF NST
! K =    NODES INDEX
! N0 =   BOUNDARY NODE TO BE ADDED TO NODES
! INDF = IADJ INDEX OF THE FIRST NEIGHBOR OF N0
!
      nn = n
!
! SEARCH FOR A BOUNDARY NODE
!
      do nst = 1 , nn
         indl = iend(nst)
         if ( iadj(indl)==0 ) goto 100
      enddo
!
! NO BOUNDARY NODE EXISTS
!
      nb = 0
      nt = 2*(nn-2)
      na = 3*(nn-2)
      return
!
! NST IS THE FIRST BOUNDARY NODE ENCOUNTERED.  INITIALIZE
!   FOR BOUNDARY TRAVERSAL.
!
 100  nodes(1) = nst
      k = 1
      n0 = nst
      do
!
! TRAVERSE THE BOUNDARY IN COUNTERCLOCKWISE ORDER
!
         indf = 1
         if ( n0>1 ) indf = iend(n0-1) + 1
         n0 = iadj(indf)
         if ( n0==nst ) then
!
! TERMINATION
!
            nb = k
            nt = 2*(nn-1) - nb
            na = 3*(nn-1) - nb
            exit
         else
            k = k + 1
            nodes(k) = n0
         endif
      enddo
      end subroutine bnodes
!*==CIRCLE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine circle(n,x,y)
!*--********************************************************************
!A INPUT  - N
!A OUTPUT - X
!A OUTPUT - Y
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   TRPLOT
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DTH      I        NM1      TH
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n
      real x(n) , y(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE COMPUTES THE COORDINATES OF A SEQUENCE
! OF EQUISPACED POINTS ON A UNIT CIRCLE.  AN (N-1)-SIDED
! POLYGONAL APPROXIMATION TO THE CIRCLE MAY BE PLOTTED BY
! CONNECTING (X(I),Y(I)) TO (X(I+1),Y(I+1)) FOR I = 1,...,
! N-1.
!
! INPUT PARAMETERS -   N - NUMBER OF POINTS.  N .GE. 2.
!
!                    X,Y - VECTORS OF LENGTH .GE. N.
!
! N IS NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - X,Y - COORDINATES OF POINTS ON THE
!                           CIRCLE WHERE (X(N),Y(N)) =
!                           (X(1),Y(1)) = (1,0).
!
! MODULES REFERENCED BY CIRCLE - NONE
!
! INTRINSIC FUNCTIONS CALLED BY CIRCLE - ATAN, FLOAT, COS,
!                                        SIN
!
!***********************************************************
!
      integer nm1 , i
      real dth , th
!
! LOCAL PARAMETERS -
!
! NM1 = N - 1
! I =   DO-LOOP INDEX
! DTH = ANGLE BETWEEN ADJACENT POINTS
! TH =  POLAR COORDINATE ANGLE
!
      nm1 = n - 1
      if ( nm1<1 ) return
      dth = 8.*atan(1.)/float(nm1)
      do i = 1 , nm1
         th = float(i-1)*dth
         x(i) = cos(th)
         y(i) = sin(th)
      enddo
      x(n) = x(1)
      y(n) = y(1)
      end subroutine circle
!*==CONSTR.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine constr(xk,yk,zk,cx,sx,cy,sy)
!*--********************************************************************
!A INPUT  - XK
!A INPUT  - YK
!A INPUT  - ZK
!A OUTPUT - CX
!A OUTPUT - SX
!A OUTPUT - CY
!A OUTPUT - SY
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADG    GRADL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  *** NONE ****
! uses PARAMs *** NONE ****
!*++********************************************************************
      real xk , yk , zk , cx , sx , cy , sy
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE CONSTRUCTS THE ELEMENTS OF A 3 BY 3
! ORTHOGONAL MATRIX R WHICH ROTATES A POINT (XK,YK,ZK) ON
! THE UNIT SPHERE TO THE NORTH POLE, I.E.
!
!      (XK)     (CY  0 -SY)   (1   0   0)   (XK)     (0)
!  R * (YK)  =  ( 0  1   0) * (0  CX -SX) * (YK)  =  (0)
!      (ZK)     (SY  0  CY)   (0  SX  CX)   (ZK)     (1)
!
! INPUT PARAMETERS - XK,YK,ZK - COMPONENTS OF A UNIT VECTOR
!                               TO BE ROTATED TO (0,0,1).
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - CX,SX,CY,SY - ELEMENTS OF R -- CX,SX
!                                   DEFINE A ROTATION ABOUT
!                                   THE X-AXIS AND CY,SY DE-
!                                   FINE A ROTATION ABOUT
!                                   THE Y-AXIS.
!
! MODULES REFERENCED BY CONSTR - NONE
!
! INTRINSIC FUNCTION CALLED BY CONSTR - SQRT
!
!***********************************************************
!
      cy = sqrt(yk*yk+zk*zk)
      sy = xk
      if ( cy==0. ) then
!
! (XK,YK,ZK) LIES ON THE X-AXIS
!
         cx = 1.
         sx = 0.
         goto 99999
      endif
      cx = zk/cy
      sx = yk/cy
      return
99999 end subroutine constr
!*==COVSPH.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine covsph(kk,node,iadj,iend)
!*--********************************************************************
!A INPUT  - KK
!A INPUT  - NODE
!A OUTPUT - IADJ
!A OUTPUT - IEND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   ADNODE
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  INDX     K        KEND     ND       NEXT
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer kk , node , iadj(*) , iend(kk)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE CONNECTS AN EXTERIOR NODE KK TO ALL
! BOUNDARY NODES OF A TRIANGULATION OF KK-1 POINTS ON THE
! UNIT SPHERE, PRODUCING A TRIANGULATION WHICH COVERS THE
! SPHERE.  IADJ AND IEND ARE UPDATED WITH THE ADDITION OF
! NODE KK, BUT NO OPTIMIZATION OF THE MESH IS PERFORMED.
! ALL BOUNDARY NODES MUST BE VISIBLE FROM KK.
!
! INPUT PARAMETERS -   KK - INDEX OF THE EXTERIOR NODE TO
!                           BE ADDED.  KK .GE. 4.
!
!                    NODE - BOUNDARY NODE INDEX IN THE
!                           RANGE 1,...,KK-1.
!
!                    IADJ - SET OF ADJACENCY LISTS FOR
!                           NODES 1,...,KK-1.
!
!                    IEND - POINTERS TO THE ENDS OF ADJA-
!                           CENCY LISTS IN IADJ FOR NODES
!                           1,...,KK-1.
!
! KK AND NODE ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  ALL NODES ARE
!                                 INTERIOR.
!
! MODULES REFERENCED BY COVSPH - NONE
!
!***********************************************************
!
      integer k , nd , next , kend , indx
!
! LOCAL PARAMETERS -
!
! K,ND = LOCAL COPIES OF KK AND NODE
! NEXT = BOUNDARY NODE TO BE CONNECTED TO K
! KEND = IADJ INDEX OF THE LAST NEIGHBOR OF K
! INDX = IADJ INDEX
!
      k = kk
      nd = node
!
! INITIALIZATION
!
      next = nd
      kend = iend(k-1)
      do
!
! WALK ALONG THE BOUNDARY CONNECTING NODE K AND NEXT.  K
!   K REPLACES 0 AS THE LAST NEIGHBOR OF NEXT.
!
         kend = kend + 1
         iadj(kend) = next
         indx = iend(next)
         iadj(indx) = k
         next = iadj(indx-1)
         if ( next==nd ) then
            iend(k) = kend
            exit
         endif
      enddo
      end subroutine covsph
!*==DELETE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine delete(n,nout1,nout2,iadj,iend,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - NOUT1
!A INPUT  - NOUT2
!A OUTPUT - IADJ
!A OUTPUT - IEND
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       INDEX    SHIFTD
! called by   AA0001   AA0003
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IMAX     IND12    IND1F    IND1L    IND21    IND2F    IND2L    INDFP2   INDLM3   INDN0    INDNF    INDNL    IO1      IO2      IOUT1    IOUT2    ITEMP    NEWBD    NF       NL       NN
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      !integer index
!*** End of declarations inserted by SPAG
      integer n , nout1 , nout2 , iadj(*) , iend(n) , ier
      !external index
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE DELETES A BOUNDARY EDGE FROM A TRIANGU-
! LATION OF A SET OF POINTS ON THE UNIT SPHERE.  IT MAY BE
! NECESSARY TO FORCE CERTAIN EDGES TO BE PRESENT BEFORE
! CALLING DELETE (SEE SUBROUTINE EDGE).  NOTE THAT SUBROU-
! TINES EDGE, TRFIND, AND THE ROUTINES WHICH CALL TRFIND
! (ADNODE, UNIF, INTRC1, AND INTRC0) SHOULD NOT BE CALLED
! FOLLOWING A DELETION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE TRIAN-
!                            GULATION.
!
!              NOUT1,NOUT2 - PAIR OF ADJACENT NODES ON THE
!                            BOUNDARY DEFINING THE ARC TO
!                            BE REMOVED.  NOUT2 MUST BE THE
!                            LAST NONZERO NEIGHBOR OF NOUT1.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                IADJ,IEND - DATA STRUCTURE DEFINING THE
!                            TRIANGULATION (SEE SUBROUTINE
!                            TRMESH).
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE REMOVAL
!                                 OF THE ARC NOUT1-NOUT2
!                                 IF IER .EQ. 0.
!
!                           IER - ERROR INDICATOR
!                                 IER = 0 IF NO ERRORS WERE
!                                         ENCOUNTERED.
!                                 IER = 1 IF NOUT1 OR NOUT2
!                                         IS NOT ON THE
!                                         BOUNDARY.
!                                 IER = 2 IF NOUT1 OR NOUT2
!                                         HAS ONLY 2 NONZERO
!                                         NEIGHBORS.
!                                 IER = 3 IF NOUT2 IS NOT
!                                         THE LAST NEIGHBOR
!                                         OF NOUT1.
!                                 IER = 4 IF A DELETION
!                                         WOULD DIVIDE THE
!                                         MESH INTO TWO
!                                         REGIONS.
!
! MODULES REFERENCED BY DELETE - SHIFTD, INDEX
!
!***********************************************************
!
      integer nn , iout1 , iout2 , io1 , io2 , ind12 , ind21 , itemp ,  &
     &        ind1f , ind1l , ind2f , ind2l , newbd , indnf , indnl ,   &
     &        indn0 , indfp2 , indlm3 , nf , nl , i , imax
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! IOUT1,IOUT2 = LOCAL COPIES OF NOUT1 AND NOUT2
! IO1,IO2 =     NOUT1,NOUT2 IN ORDER OF INCREASING MAGNITUDE
! IND12 =       INDEX OF IO2 IN THE ADJACENCY LIST FOR IO1
! IND21 =       INDEX OF IO1 IN THE ADJACENCY LIST FOR IO2
! ITEMP =       TEMPORARY STORAGE LOCATION FOR PERMUTATIONS
! IND1F =       IADJ INDEX OF THE FIRST NEIGHBOR OF IO1
! IND1L =       IADJ INDEX OF THE LAST NEIGHBOR OF IO1
! IND2F =       IADJ INDEX OF THE FIRST NEIGHBOR OF IO2
! IND2L =       IADJ INDEX OF THE LAST NEIGHBOR OF IO2
! NEWBD =       THE NEIGHBOR COMMON TO NOUT1 AND NOUT2
! INDNF =       IADJ INDEX OF THE FIRST NEIGHBOR OF NEWBD
! INDNL =       IADJ INDEX OF THE LAST NEIGHBOR OF NEWBD
! INDN0 =       INDEX OF 0 IN THE ADJACENCY LIST FOR NEWBD
!                 BEFORE PERMUTING THE NEIGHBORS
! INDFP2 =      INDNF + 2
! INDLM3 =      INDNL - 3
! NF,NL =       BOUNDS ON THE PORTION OF IADJ TO BE SHIFTED
! I =           DO-LOOP INDEX
! IMAX =        UPPER BOUND ON DO-LOOP FOR SHIFTING IEND
!
      nn = n
      iout1 = nout1
      iout2 = nout2
!
! INITIALIZE INDICES
!
      ind1f = 1
      if ( iout1>1 ) ind1f = iend(iout1-1) + 1
      ind1l = iend(iout1)
      ind2f = 1
      if ( iout2>1 ) ind2f = iend(iout2-1) + 1
      ind2l = iend(iout2)
      newbd = iadj(ind1l-2)
      indn0 = indexr(newbd,iout2,iadj,iend)
      indnl = iend(newbd)
!
! ORDER VERTICES SUCH THAT THE NEIGHBORS OF IO1 PRECEDE
!   THOSE OF IO2
!
      if ( iout1>iout2 ) then
         io1 = iout2
         io2 = iout1
         ind12 = ind2f
         ind21 = ind1l - 1
      else
         io1 = iout1
         io2 = iout2
         ind12 = ind1l - 1
         ind21 = ind2f
      endif
!
! CHECK FOR ERRORS
!
      if ( (iadj(ind1l)/=0) .or. (iadj(ind2l)/=0) ) then
!
! ONE OF THE VERTICES IS NOT ON THE BOUNDARY
!
         ier = 1
         return
      elseif ( (ind1l-ind1f<=2) .or. (ind2l-ind2f<=2) ) then
!
! ONE OF THE VERTICES HAS ONLY TWO NONZERO NEIGHBORS -- THE
!   TRIANGULATION WOULD BE DESTROYED BY A DELETION
!
         ier = 2
         return
      elseif ( iadj(ind1l-1)/=iout2 ) then
!
! NOUT2 IS NOT THE LAST NONZERO NEIGHBOR OF NOUT1
!
         ier = 3
         return
      else
         if ( iadj(indnl)==0 ) then
!
! A DELETION WOULD DIVIDE THE MESH INTO TWO REGIONS
!   CONNECTED AT A SINGLE NODE
!
            ier = 4
            goto 99999
         else
!
! DELETE THE EDGE IO1-IO2 AND MAKE NEWBD A BOUNDARY NODE
!
            if ( newbd<io1 ) then
!
! THE VERTICES ARE ORDERED NEWBD, IO1, IO2.
! DELETE IO2 AS A NEIGHBOR OF IO1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD.
!
               indn0 = indn0 + 1
               nf = indn0
               nl = ind12 - 1
               if ( nf<=nl ) call shiftd(nf,nl,1,iadj)
               imax = io1 - 1
               do i = newbd , imax
                  iend(i) = iend(i) + 1
               enddo
            elseif ( newbd<io2 ) then
!
! THE VERTICES ARE ORDERED IO1, NEWBD, IO2.
! DELETE IO2 AS A NEIGHBOR OF IO1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD.
!
               nf = ind12 + 1
               nl = indn0
               call shiftd(nf,nl,-1,iadj)
               imax = newbd - 1
               do i = io1 , imax
                  iend(i) = iend(i) - 1
               enddo
            else
!
! THE VERTICES ARE ORDERED IO1, IO2, NEWBD.
! DELETE IO2 AS A NEIGHBOR OF IO1.
!
               nf = ind12 + 1
               nl = ind21 - 1
               call shiftd(nf,nl,-1,iadj)
               imax = io2 - 1
               do i = io1 , imax
                  iend(i) = iend(i) - 1
               enddo
!
! DELETE IO1 AS A NEIGHBOR OF IO2
!
               nf = nl + 2
               nl = indn0
               call shiftd(nf,nl,-2,iadj)
               imax = newbd - 1
               do i = io2 , imax
                  iend(i) = iend(i) - 2
               enddo
!
! SHIFT THE BOTTOM OF IADJ UP 1 LEAVING ROOM FOR 0 AS A
!   NEIGHBOR OF NEWBD
!
               indn0 = indn0 - 1
               nf = nl + 1
               nl = iend(nn)
               if ( nf<=nl ) call shiftd(nf,nl,-1,iadj)
               do i = newbd , nn
                  iend(i) = iend(i) - 1
               enddo
               goto 50
            endif
!
! DELETE IO1 AS A NEIGHBOR OF IO2
!
            nf = ind21 + 1
            nl = iend(nn)
            call shiftd(nf,nl,-1,iadj)
            do i = io2 , nn
               iend(i) = iend(i) - 1
            enddo
         endif
!
! PERMUTE THE NEIGHBORS OF NEWBD WITH END-AROUND SHIFTS SO
!   THAT 0 IS THE LAST NEIGHBOR
!
 50      indnf = 1
         if ( newbd>1 ) indnf = iend(newbd-1) + 1
         indnl = iend(newbd)
         if ( indn0-indnf>=indnl-indn0 ) then
!
! SHIFT DOWNWARD
!
            if ( indn0/=indnl ) then
               if ( indn0<indnl-1 ) then
                  indlm3 = indnl - 3
                  if ( indn0<=indlm3 ) then
                     do i = indn0 , indlm3
                        itemp = iadj(indnl)
                        call shiftd(indnf,indnl-1,1,iadj)
                        iadj(indnf) = itemp
                     enddo
                  endif
!
! THE LAST SHIFT IS BY 2
!
                  itemp = iadj(indnl-1)
                  call shiftd(indnf,indlm3,2,iadj)
                  iadj(indnf+1) = iadj(indnl)
                  iadj(indnf) = itemp
               else
                  call shiftd(indnf,indnl-2,1,iadj)
                  iadj(indnf) = iadj(indnl)
               endif
            endif
!
! SHIFT UPWARD
!
         elseif ( indn0>indnf ) then
            indfp2 = indnf + 2
            if ( indn0>=indfp2 ) then
               do i = indfp2 , indn0
                  itemp = iadj(indnf)
                  call shiftd(indnf+1,indnl,-1,iadj)
                  iadj(indnl) = itemp
               enddo
            endif
!
! THE LAST SHIFT IS BY 2
!
            itemp = iadj(indnf)
            call shiftd(indfp2,indnl,-2,iadj)
            iadj(indnl-1) = itemp
         else
            call shiftd(indnf+1,indnl,-1,iadj)
         endif
      endif
!
! INSERT 0 AS THE LAST NEIGHBOR OF NEWBD
!
      iadj(indnl) = 0
      ier = 0
      return
99999 end subroutine delete
!*==EDGE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine edge(in1,in2,x,y,z,lwk,iwk,iadj,iend,ier)
!*--********************************************************************
!A INPUT  - IN1
!A INPUT  - IN2
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A OUTPUT - LWK
!A OUTPUT - IWK
!A INPUT  - IADJ
!A INPUT  - IEND
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       SWAP     SWPTST
! called by   AA0001   AA0003
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        INDF     INDL     INDX     IO1      IO2      IWC      IWCM1    IWCP1    IWEND    IWF      IWL      LFT      N0       N1       N1LST    N2       NEXT     NL       NR       SWP      X0       X1       X2       XA       XB       XP
!             Y0       Y1       Y2       YA       YB       YP       Z0       Z1       Z2       ZA       ZB       ZP
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real xa , xb , xp , ya , yb , yp , za , zb , zp
!*** End of declarations inserted by SPAG
      !logical swptst
      integer in1 , in2 , lwk , iwk(2,*) , iadj(*) , iend(*) , ier
      real x(*) , y(*) , z(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N NODES AND A PAIR OF NODAL
! INDICES IN1 AND IN2, THIS ROUTINE SWAPS ARCS AS NECESSARY
! TO FORCE IN1 AND IN2 TO BE ADJACENT.  ONLY ARCS WHICH
! INTERSECT IN1-IN2 ARE SWAPPED OUT.  IF A THIESSEN TRIANGU-
! LATION IS INPUT, THE RESULTING TRIANGULATION IS AS CLOSE
! AS POSSIBLE TO A THIESSEN TRIANGULATION IN THE SENSE THAT
! ALL ARCS OTHER THAN IN1-IN2 ARE LOCALLY OPTIMAL.
!   A SEQUENCE OF CALLS TO EDGE MAY BE USED TO FORCE THE
! PRESENCE OF A SET OF EDGES DEFINING THE BOUNDARY OF A NON-
! CONVEX REGION.  SUBSEQUENT DELETION OF EDGES OUTSIDE THIS
! REGION (BY SUBROUTINE DELETE) RESULTS IN A NONCONVEX TRI-
! ANGULATION WHICH MAY SERVE AS A FINITE ELEMENT GRID.
! (EDGE SHOULD NOT BE CALLED AFTER A CALL TO DELETE.)  IF,
! ON THE OTHER HAND, INTERPOLATION IS TO BE PERFORMED IN THE
! NONCONVEX REGION, EDGES MUST NOT BE DELETED, BUT IT IS
! STILL ADVANTAGEOUS TO HAVE THE NONCONVEX BOUNDARY PRESENT
! IF IT IS DESIRABLE THAT INTERPOLATED VALUES BE INFLUENCED
! BY THE GEOMETRY.  NOTE THAT SUBROUTINE GETNP WHICH IS USED
! TO SELECT THE NODES ENTERING INTO LOCAL DERIVATIVE ESTI-
! MATES WILL NOT NECESSARILY RETURN CLOSEST NODES IF THE
! TRIANGULATION HAS BEEN RENDERED NONOPTIMAL BY A CALL TO
! EDGE.  HOWEVER, THE EFFECT WILL BE MERELY TO FURTHER EN-
! HANCE THE INFLUENCE OF THE NONCONVEX GEOMETRY ON INTERPO-
! LATED VALUES.
!
! INPUT PARAMETERS - IN1,IN2 - INDICES (OF X,Y,Z) IN THE
!                              RANGE 1,...,N DEFINING A PAIR
!                              OF NODES TO BE CONNECTED BY
!                              AN ARC.
!
!                      X,Y,Z - N-VECTORS CONTAINING CARTE-
!                              SIAN COORDINATES OF THE
!                              NODES.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                        LWK - NUMBER OF COLUMNS RESERVED
!                              FOR IWK.  THIS MUST BE AT
!                              LEAST NI -- THE NUMBER OF
!                              ARCS WHICH INTERSECT IN1-IN2.
!                              (NI IS BOUNDED BY N-3).
!
!                        IWK - INTEGER WORK ARRAY DIMENSION-
!                              ED 2 BY LWK (OR VECTOR OF
!                              LENGTH .GE. 2*LWK).
!
!                  IADJ,IEND - DATA STRUCTURE DEFINING THE
!                              TRIANGULATION.  SEE SUBROU-
!                              TINE TRMESH.
!
! OUTPUT PARAMETERS - LWK - NUMBER OF IWK COLUMNS REQUIRED
!                           IF IER = 0 OR IER = 2.  LWK = 0
!                           IFF IN1 AND IN2 WERE ADJACENT
!                           ON INPUT.
!
!                     IWK - CONTAINS THE INDICES OF THE END-
!                           POINTS OF THE NEW ARCS OTHER
!                           THAN IN1-IN2 UNLESS IER .GT. 0
!                           OR LWK = 0.  NEW ARCS TO THE
!                           LEFT OF IN1->IN2 ARE STORED IN
!                           THE FIRST K-1 COLUMNS (LEFT POR-
!                           TION OF IWK), COLUMN K CONTAINS
!                           ZEROS, AND NEW ARCS TO THE RIGHT
!                           OF IN1->IN2 OCCUPY COLUMNS K+1,
!                           ...,LWK.  (K CAN BE DETERMINED
!                           BY SEARCHING IWK FOR THE ZEROS.)
!
!               IADJ,IEND - UPDATED IF NECESSARY TO REFLECT
!                           THE PRESENCE OF AN ARC CONNECT-
!                           ING IN1 AND IN2, UNALTERED IF
!                           IER .NE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE EN-
!                                   COUNTERED.
!                           IER = 1 IF IN1 .LT. 1, IN2 .LT.
!                                   1, IN1 = IN2, OR LWK
!                                   .LT. 0 ON INPUT.
!                           IER = 2 IF MORE SPACE IS REQUIR-
!                                   ED IN IWK.  SEE LWK.
!                           IER = 3 IF IN1 AND IN2 COULD NOT
!                                   BE CONNECTED DUE TO AN
!                                   INVALID DATA STRUCTURE.
!
! MODULES REFERENCED BY EDGE - SWAP, INDEX, SHIFTD, SWPTST
!
!***********************************************************
!
      integer n1 , n2 , iwend , iwl , indf , indx , n1lst , nl , nr ,   &
     &        next , iwf , lft , n0 , iwc , iwcp1 , iwcm1 , i , io1 ,   &
     &        io2 , indl
      real x1 , y1 , z1 , x2 , y2 , z2 , x0 , y0 , z0
      logical swp , left
!
! LOCAL PARAMETERS -
!
! N1,N2 =    LOCAL COPIES OF IN1 AND IN2 OR NODES OPPOSITE
!              AN ARC IO1-IO2 TO BE TESTED FOR A SWAP IN
!              THE OPTIMIZATION LOOPS
! IWEND =    INPUT OR OUTPUT VALUE OF LWK
! IWL =      IWK (COLUMN) INDEX OF THE LAST (RIGHTMOST) ARC
!              WHICH INTERSECTS IN1->IN2
! INDF =     IADJ INDEX OF THE FIRST NEIGHBOR OF IN1 OR IO1
! INDX =     IADJ INDEX OF A NEIGHBOR OF IN1, NL, OR IO1
! N1LST =    LAST NEIGHBOR OF IN1
! NL,NR =    ENDPOINTS OF AN ARC WHICH INTERSECTS IN1-IN2
!              WITH NL LEFT IN1->IN2
! NEXT =     NODE OPPOSITE NL->NR
! IWF =      IWK (COLUMN) INDEX OF THE FIRST (LEFTMOST) ARC
!              WHICH INTERSECTS IN1->IN2
! LFT =      FLAG USED TO DETERMINE IF A SWAP RESULTS IN THE
!              NEW ARC INTERSECTING IN1-IN2 -- LFT = 0 IFF
!              N0 = IN1, LFT = -1 IMPLIES N0 LEFT IN1->IN2,
!              AND LFT = 1 IMPLIES N0 LEFT IN2->IN1
! N0 =       NODE OPPOSITE NR->NL
! IWC =      IWK INDEX BETWEEN IWF AND IWL -- NL->NR IS
!              STORED IN IWK(1,IWC)->IWK(2,IWC)
! IWCP1 =    IWC + 1
! IWCM1 =    IWC - 1
! I =        DO-LOOP INDEX AND COLUMN INDEX FOR IWK
! IO1,IO2 =  ENDPOINTS OF AN ARC TO BE TESTED FOR A SWAP IN
!              THE OPTIMIZATION LOOPS
! INDL =     IADJ INDEX OF THE LAST NEIGHBOR OF IO1
! X1,Y1,Z1 = COORDINATES OF IN1
! X2,Y2,Z2 = COORDINATES OF IN2
! X0,Y0,Z0 = COORDINATES OF N0
! SWP =      FLAG SET TO .TRUE. IFF A SWAP OCCURS IN AN OPT-
!              IMIZATION LOOP
! LEFT =     STATEMENT FUNCTION WHICH RETURNS THE VALUE
!              .TRUE. IFF (XP,YP,ZP) IS ON OR TO THE LEFT OF
!              THE VECTOR (XA,YA,ZA)->(XB,YB,ZB)
!
      left(xa,ya,za,xb,yb,zb,xp,yp,zp) = xp*(ya*zb-yb*za)               &
     &   - yp*(xa*zb-xb*za) + zp*(xa*yb-xb*ya)>=0.
!
! STORE IN1, IN2, AND LWK IN LOCAL VARIABLES AND CHECK FOR
!   ERRORS.
!
      n1 = in1
      n2 = in2
      iwend = lwk
      if ( n1<1 .or. n2<1 .or. n1==n2 .or. iwend<0 ) then
!
! PARAMETER OUT OF RANGE
!
         ier = 1
         return
      else
!
! STORE THE COORDINATES OF N1 AND N2 AND INITIALIZE IWL.
!
         x1 = x(n1)
         y1 = y(n1)
         z1 = z(n1)
         x2 = x(n2)
         y2 = y(n2)
         z2 = z(n2)
         iwl = 0
!
! SET NR AND NL TO ADJACENT NEIGHBORS OF N1 SUCH THAT
!   NR LEFT N2->N1 AND NL LEFT N1->N2.
!
!   SET INDF AND INDX TO THE INDICES OF THE FIRST AND LAST
!     NEIGHBORS OF N1 AND SET N1LST TO THE LAST NEIGHBOR.
!
         indf = 1
         if ( n1>1 ) indf = iend(n1-1) + 1
         indx = iend(n1)
         n1lst = iadj(indx)
         if ( n1lst==0 ) indx = indx - 1
         if ( n1lst/=0 ) then
!
!   N1 IS AN INTERIOR NODE.  LOOP THROUGH THE NEIGHBORS NL
!     IN REVERSE ORDER UNTIL NL LEFT N1->N2.
!
            nl = n1lst
            do while ( .not.(left(x1,y1,z1,x2,y2,z2,x(nl),y(nl),z(nl))) &
     &                 )
               indx = indx - 1
               nl = iadj(indx)
               if ( indx<=indf ) then
!
!   NL IS THE FIRST NEIGHBOR OF N1.  SET NR TO THE LAST
!     NEIGHBOR AND TEST FOR AN ARC N1-N2.
!
                  nr = n1lst
                  if ( nl/=n2 ) goto 100
                  goto 700
               endif
            enddo
         endif
         do
!
!   NL = IADJ(INDX) LEFT N1->N2 AND INDX .GT. INDF.  SET
!     NR TO THE PRECEDING NEIGHBOR OF N1.
!
            indx = indx - 1
            nr = iadj(indx)
            if ( left(x2,y2,z2,x1,y1,z1,x(nr),y(nr),z(nr)) ) then
!
!   SET NL TO THE NEIGHBOR FOLLOWING NR AND TEST FOR AN ARC
!     N1-N2.
!
               nl = iadj(indx+1)
               if ( nl/=n2 .and. nr/=n2 ) exit
               goto 700
            elseif ( indx<=indf ) then
!
!   SET NL AND NR TO THE FIRST AND LAST NEIGHBORS OF N1 AND
!     TEST FOR AN INVALID DATA STRUCTURE (N1 CANNOT BE A
!     BOUNDARY NODE AND CANNOT BE ADJACENT TO N2).
!
               nl = nr
               nr = n1lst
               if ( nr/=0 .and. nr/=n2 ) exit
!
! INVALID TRIANGULATION DATA STRUCTURE
!
               ier = 3
               goto 99999
            endif
         enddo
      endif
!
! STORE THE ORDERED SEQUENCE OF INTERSECTING EDGES NL->NR IN
!   IWK(1,IWL)->IWK(2,IWL).
!
 100  iwl = iwl + 1
      if ( iwl<=iwend ) iwk(1,iwl) = nl
      if ( iwl<=iwend ) iwk(2,iwl) = nr
!
!   SET NEXT TO THE NEIGHBOR OF NL WHICH FOLLOWS NR.
!
      indx = iend(nl)
      if ( iadj(indx)/=nr ) then
         do
!
!   NR IS NOT THE LAST NEIGHBOR OF NL.  LOOP THROUGH THE
!     NEIGHBORS IN REVERSE ORDER.
!
            indx = indx - 1
            if ( iadj(indx)==nr ) exit
         enddo
      else
!
!   NR IS THE LAST NEIGHBOR OF NL.  SET NEXT TO THE FIRST
!     NEIGHBOR.
!
         indx = 0
         if ( nl/=1 ) indx = iend(nl-1)
      endif
!
!   STORE NEXT, TEST FOR AN INVALID TRIANGULATION (NL->NR
!     CANNOT BE A BOUNDARY EDGE), AND TEST FOR TERMINATION
!     OF THE LOOP.
!
      next = iadj(indx+1)
      if ( next==0 ) then
         ier = 3
         goto 99999
      elseif ( next==n2 ) then
!
! IWL IS THE NUMBER OF ARCS WHICH INTERSECT N1-N2.  STORE
!   LWK AND TEST FOR SUFFICIENT SPACE.
!
         lwk = iwl
         if ( iwl>iwend ) then
!
! INSUFFICIENT SPACE IN IWK
!
            ier = 2
            return
         else
            iwend = iwl
!
! INITIALIZE FOR EDGE SWAPPING LOOP -- ALL POSSIBLE SWAPS
!   ARE APPLIED (EVEN IF THE NEW ARC AGAIN INTERSECTS
!   N1-N2), ARCS TO THE LEFT OF N1->N2 ARE STORED IN THE
!   LEFT PORTION OF IWK, AND ARCS TO THE RIGHT ARE STORED IN
!   THE RIGHT PORTION.  IWF AND IWL INDEX THE FIRST AND LAST
!   INTERSECTING ARCS.
!
            ier = 0
            iwf = 1
         endif
      else
!
!   SET NL OR NR TO NEXT.
!
         if ( left(x1,y1,z1,x2,y2,z2,x(next),y(next),z(next)) ) then
            nl = next
         else
            nr = next
         endif
         goto 100
      endif
!
! TOP OF LOOP -- SET N0 TO N1 AND NL->NR TO THE FIRST EDGE.
!   IWC POINTS TO THE ARC CURRENTLY BEING PROCESSED.  LFT
!   .LE. 0 IFF N0 LEFT N1->N2.
!
 200  lft = 0
      n0 = n1
      x0 = x1
      y0 = y1
      z0 = z1
      nl = iwk(1,iwf)
      nr = iwk(2,iwf)
      iwc = iwf
!
!   SET NEXT TO THE NODE OPPOSITE NL->NR UNLESS IWC IS THE
!     LAST ARC.
!
 300  do while ( iwc/=iwl )
         iwcp1 = iwc + 1
         next = iwk(1,iwcp1)
         if ( next/=nl ) then
!
!   NEXT LEFT N1->N2, NEXT .NE. N2, AND IWC .LT. IWL.
!     TEST FOR A POSSIBLE SWAP.
!
            if ( .not.left(x(nl),y(nl),z(nl),x0,y0,z0,x(next),y(next),  &
     &           z(next)) ) goto 400
            if ( lft<=0 ) then
!
!   SWAP NL-NR FOR N0-NEXT, SHIFT COLUMNS IWF,...,IWC-1 TO
!     THE RIGHT, AND STORE N0-NEXT IN THE LEFT PORTION OF
!     IWK.
!
               call swap(next,n0,nl,nr,iadj,iend)
               i = iwc
               do while ( i/=iwf )
                  iwk(1,i) = iwk(1,i-1)
                  iwk(2,i) = iwk(2,i-1)
                  i = i - 1
               enddo
               iwk(1,iwf) = n0
               iwk(2,iwf) = next
               iwf = iwf + 1
               goto 500
            else
               if ( .not.left(x0,y0,z0,x(nr),y(nr),z(nr),x(next),y(next)&
     &              ,z(next)) ) goto 400
!
!   REPLACE NL->NR WITH NEXT->N0.
!
               call swap(next,n0,nl,nr,iadj,iend)
               iwk(1,iwc) = next
               iwk(2,iwc) = n0
               goto 500
            endif
         else
            next = iwk(2,iwcp1)
!
!   NEXT RIGHT N1->N2 AND IWC .LT. IWL.  TEST FOR A POSSIBLE
!     SWAP.
!
            if ( left(x0,y0,z0,x(nr),y(nr),z(nr),x(next),y(next),z(next)&
     &           ) ) then
               if ( lft>=0 ) then
!
!   SWAP NL-NR FOR N0-NEXT, SHIFT COLUMNS IWC+1,...,IWL TO
!     THE LEFT, AND STORE N0-NEXT IN THE RIGHT PORTION OF
!     IWK.
!
                  call swap(next,n0,nl,nr,iadj,iend)
                  do i = iwcp1 , iwl
                     iwk(1,i-1) = iwk(1,i)
                     iwk(2,i-1) = iwk(2,i)
                  enddo
                  iwk(1,iwl) = n0
                  iwk(2,iwl) = next
                  iwl = iwl - 1
                  nr = next
                  cycle
               elseif ( left(x(nl),y(nl),z(nl),x0,y0,z0,x(next),y(next),&
     &                  z(next)) ) then
!
!   REPLACE NL->NR WITH N0->NEXT.
!
                  call swap(next,n0,nl,nr,iadj,iend)
                  iwk(1,iwc) = n0
                  iwk(2,iwc) = next
                  goto 350
               endif
            endif
!
!   A SWAP IS NOT POSSIBLE.  SET N0 TO NR.
!
            n0 = nr
            x0 = x(n0)
            y0 = y(n0)
            z0 = z(n0)
            lft = 1
         endif
!
!   ADVANCE TO THE NEXT ARC.
!
 350     nr = next
         iwc = iwc + 1
      enddo
!
!   N2 IS OPPOSITE NL->NR (IWC = IWL).
!
      if ( n0==n1 ) then
!
! IWF = IWC = IWL.  SWAP OUT THE LAST ARC FOR N1-N2 AND
!   STORE ZEROS IN IWK.
!
         call swap(n2,n1,nl,nr,iadj,iend)
         iwk(1,iwc) = 0
         iwk(2,iwc) = 0
         if ( iwc==1 ) goto 600
!
! OPTIMIZATION LOOPS -- OPTIMIZE THE SET OF NEW ARCS TO THE
!   LEFT OF IN1->IN2.  THE LOOP IS REPEATED UNTIL NO SWAPS
!   ARE PERFORMED.
!
         iwcm1 = iwc - 1
         do
            swp = .false.
            do i = 1 , iwcm1
               io1 = iwk(1,i)
               io2 = iwk(2,i)
!
!   SET N1 TO THE NEIGHBOR OF IO1 WHICH FOLLOWS IO2 AND SET
!     N2 TO THE NEIGHBOR OF IO1 WHICH PRECEDES IO2.
!
               indf = 1
               if ( io1>1 ) indf = iend(io1-1) + 1
               indl = iend(io1)
               indx = indl
               if ( iadj(indx)/=io2 ) then
                  do
!
!   IO2 IS NOT THE LAST NEIGHBOR OF IO1.  LOOP THROUGH THE
!     NEIGHBORS IN REVERSE ORDER.
!
                     indx = indx - 1
                     if ( iadj(indx)==io2 ) then
                        n1 = iadj(indx+1)
                        if ( indx/=indf ) n2 = iadj(indx-1)
                        if ( indx==indf ) n2 = iadj(indl)
                        exit
                     endif
                  enddo
               else
!
!   IO2 IS THE LAST NEIGHBOR OF IO1.
!
                  n1 = iadj(indf)
                  n2 = iadj(indx-1)
               endif
!
!   TEST IO1-IO2 FOR A SWAP.
!
               if ( swptst(n1,n2,io1,io2,x,y,z) ) then
                  swp = .true.
                  call swap(n1,n2,io1,io2,iadj,iend)
                  iwk(1,i) = n1
                  iwk(2,i) = n2
               endif
            enddo
            if ( .not.(swp) ) goto 600
         enddo
      elseif ( lft<0 ) then
!
!   N0 LEFT N1->N2.  TEST FOR A POSSIBLE SWAP.
!
         if ( .not.left(x(nl),y(nl),z(nl),x0,y0,z0,x2,y2,z2) ) goto 200
!
!   SWAP NL-NR FOR N0-N2, SHIFT COLUMNS IWF,...,IWL-1 TO THE
!     RIGHT, AND STORE N0-N2 IN THE LEFT PORTION OF IWK.
!
         call swap(n2,n0,nl,nr,iadj,iend)
         i = iwl
         do
            iwk(1,i) = iwk(1,i-1)
            iwk(2,i) = iwk(2,i-1)
            i = i - 1
            if ( i<=iwf ) then
               iwk(1,iwf) = n0
               iwk(2,iwf) = n2
               iwf = iwf + 1
               goto 200
            endif
         enddo
      else
!
!   N0 RIGHT N1->N2.  TEST FOR A POSSIBLE SWAP.
!
         if ( left(x0,y0,z0,x(nr),y(nr),z(nr),x2,y2,z2) ) then
!
!   SWAP NL-NR FOR N0-N2 AND STORE N0-N2 IN THE RIGHT
!     PORTION OF IWK.
!
            call swap(n2,n0,nl,nr,iadj,iend)
            iwk(1,iwl) = n0
            iwk(2,iwl) = n2
            iwl = iwl - 1
         endif
         goto 200
      endif
!
!   A SWAP IS NOT POSSIBLE.  SET N0 TO NL.
!
 400  n0 = nl
      x0 = x(n0)
      y0 = y(n0)
      z0 = z(n0)
      lft = -1
!
!   ADVANCE TO THE NEXT ARC.
!
 500  nl = next
      iwc = iwc + 1
      goto 300
!
! TEST FOR TERMINATION.
!
 600  if ( iwc==iwend ) return
      iwcp1 = iwc + 1
      do
!
! OPTIMIZE THE SET OF NEW ARCS TO THE RIGHT OF IN1->IN2.
!
         swp = .false.
         do i = iwcp1 , iwend
            io1 = iwk(1,i)
            io2 = iwk(2,i)
!
!   SET N1 AND N2 TO THE NODES OPPOSITE IO1->IO2 AND
!     IO2->IO1, RESPECTIVELY.
!
            indf = 1
            if ( io1>1 ) indf = iend(io1-1) + 1
            indl = iend(io1)
            indx = indl
            if ( iadj(indx)/=io2 ) then
               do
!
                  indx = indx - 1
                  if ( iadj(indx)==io2 ) then
                     n1 = iadj(indx+1)
                     if ( indx/=indf ) n2 = iadj(indx-1)
                     if ( indx==indf ) n2 = iadj(indl)
                     exit
                  endif
               enddo
            else
!
               n1 = iadj(indf)
               n2 = iadj(indx-1)
            endif
!
            if ( swptst(n1,n2,io1,io2,x,y,z) ) then
               swp = .true.
               call swap(n1,n2,io1,io2,iadj,iend)
               iwk(1,i) = n1
               iwk(2,i) = n2
            endif
         enddo
         if ( .not.(swp) ) return
      enddo
!
! IN1 AND IN2 WERE ADJACENT ON INPUT.
!
 700  ier = 0
      lwk = 0
      return
99999 end subroutine edge
!*==GETNP.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine getnp(x,y,z,iadj,iend,l,npts,df,ier)
!*--********************************************************************
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - IADJ
!A OUTPUT - IEND
!A INPUT  - L
!A OUTPUT - NPTS
!A OUTPUT - DF
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0001   AA0003   GRADL    GRADL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DNB      DNP      I        INDF     INDL     INDX     LM1      N1       NB       NI       NP       X1       Y1       Z1
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer iadj(*) , iend(*) , l , npts(l) , ier
      real x(*) , y(*) , z(*) , df
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF N NODES ON THE UNIT
! SPHERE AND AN ARRAY NPTS CONTAINING THE INDICES OF L-1
! NODES ORDERED BY ANGULAR DISTANCE FROM NPTS(1), THIS SUB-
! ROUTINE SETS NPTS(L) TO THE INDEX OF THE NEXT NODE IN THE
! SEQUENCE -- THE NODE, OTHER THAN NPTS(1),...,NPTS(L-1),
! WHICH IS CLOSEST TO NPTS(1).  THUS, THE ORDERED SEQUENCE
! OF K CLOSEST NODES TO N1 (INCLUDING N1) MAY BE DETERMINED
! BY K-1 CALLS TO GETNP WITH NPTS(1) = N1 AND L = 2,3,...,K
! FOR K .GE. 2.
!   THE ALGORITHM USES THE FACT THAT, IN A THIESSEN TRIAN-
! GULATION, THE K-TH CLOSEST NODE TO A GIVEN NODE N1 IS A
! NEIGHBOR OF ONE OF THE K-1 CLOSEST NODES TO N1.
!
! INPUT PARAMETERS - X,Y,Z - VECTORS OF LENGTH N CONTAINING
!                            THE CARTESIAN COORDINATES OF
!                            THE NODES.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE TRIANGULATION.
!
!                     IEND - POINTERS TO THE ENDS OF ADJA-
!                            CENCY LISTS FOR EACH NODE IN
!                            THE TRIANGULATION.
!
!                        L - NUMBER OF NODES IN THE SEQUENCE
!                            ON OUTPUT.  2 .LE. L .LE. N.
!
!                     NPTS - ARRAY OF LENGTH .GE. L CONTAIN-
!                            ING THE INDICES OF THE L-1
!                            CLOSEST NODES TO NPTS(1) IN THE
!                            FIRST L-1 LOCATIONS.
!
! IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS OTHER THAN NPTS ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - NPTS - UPDATED WITH THE INDEX OF THE
!                            L-TH CLOSEST NODE TO NPTS(1) IN
!                            POSITION L UNLESS IER = 1.
!
!                       DF - INCREASING FUNCTION (NEGATIVE
!                            COSINE) OF THE ANGULAR DISTANCE
!                            BETWEEN NPTS(1) AND NPTS(L)
!                            UNLESS IER = 1.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE EN-
!                                    COUNTERED.
!                            IER = 1 IF L IS OUT OF RANGE.
!
! MODULES REFERENCED BY GETNP - NONE
!
! INTRINSIC FUNCTION CALLED BY GETNP - IABS
!
!***********************************************************
!
      integer lm1 , n1 , i , ni , np , indf , indl , indx , nb
      real x1 , y1 , z1 , dnp , dnb
!
! LOCAL PARAMETERS -
!
! LM1 =      L - 1
! N1 =       NPTS(1)
! I =        NPTS INDEX AND DO-LOOP INDEX
! NI =       NPTS(I)
! NP =       CANDIDATE FOR NPTS(L)
! INDF =     IADJ INDEX OF THE FIRST NEIGHBOR OF NI
! INDL =     IADJ INDEX OF THE LAST NEIGHBOR OF NI
! INDX =     IADJ INDEX IN THE RANGE INDF,...,INDL
! NB =       NEIGHBOR OF NI AND CANDIDATE FOR NP
! X1,Y1,Z1 = COORDINATES OF N1
! DNP,DNB =  NEGATIVE COSINES OF THE ANGULAR DISTANCES FROM
!              N1 TO NP AND TO NB, RESPECTIVELY
!
      lm1 = l - 1
      if ( lm1<1 ) then
!
! L IS OUT OF RANGE
!
         ier = 1
         goto 99999
      endif
      ier = 0
      n1 = npts(1)
      x1 = x(n1)
      y1 = y(n1)
      z1 = z(n1)
!
! MARK THE ELEMENTS OF NPTS
!
      do i = 1 , lm1
         ni = npts(i)
         iend(ni) = -iend(ni)
      enddo
!
! CANDIDATES FOR NP = NPTS(L) ARE THE UNMARKED NEIGHBORS
!   OF NODES IN NPTS.  DNP IS INITIALIZED TO -COS(PI) --
!   THE MAXIMUM DISTANCE.
!
      dnp = 1.
!
! LOOP ON NODES NI IN NPTS
!
      do i = 1 , lm1
         ni = npts(i)
         indf = 1
         if ( ni>1 ) indf = iabs(iend(ni-1)) + 1
         indl = -iend(ni)
!
! LOOP ON NEIGHBORS NB OF NI
!
         do indx = indf , indl
            nb = iadj(indx)
            if ( nb/=0 .and. iend(nb)>=0 ) then
!
! NB IS AN UNMARKED NEIGHBOR OF NI.  REPLACE NP IF NB IS
!   CLOSER TO N1.
!
               dnb = -(x(nb)*x1+y(nb)*y1+z(nb)*z1)
               if ( dnb<dnp ) then
                  np = nb
                  dnp = dnb
               endif
            endif
         enddo
      enddo
      npts(l) = np
      df = dnp
!
! UNMARK THE ELEMENTS OF NPTS
!
      do i = 1 , lm1
         ni = npts(i)
         iend(ni) = -iend(ni)
      enddo
      return
99999 end subroutine getnp
!*==GIVENS.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine givens(a,b,c,s)
!*--********************************************************************
!A OUTPUT - A
!A OUTPUT - B
!A OUTPUT - C
!A OUTPUT - S
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADL    GRADL    QSHEP2   QSHEP3
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  AA       BB       R        U        V
! uses PARAMs *** NONE ****
!*++********************************************************************
      real a , b , c , s
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --
!     ( C  S)
! G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND
!     (-S  C)
! ENTRY OF THE 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS
! IS NORMALLY FOLLOWED BY A CALL TO ROTATE WHICH APPLIES
! THE TRANSFORMATION TO A 2 BY N MATRIX.  THIS ROUTINE WAS
! TAKEN FROM LINPACK.
!
! INPUT PARAMETERS - A,B - COMPONENTS OF THE 2-VECTOR TO BE
!                          ROTATED.
!
! OUTPUT PARAMETERS -  A - OVERWRITTEN BY R = +/-SQRT(A*A
!                          + B*B)
!
!                      B - OVERWRITTEN BY A VALUE Z WHICH
!                          ALLOWS C AND S TO BE RECOVERED
!                          AS FOLLOWS -
!                          C = SQRT(1-Z*Z), S=Z IF ABS(Z)
!                              .LE. 1.
!                          C = 1/Z, S = SQRT(1-C*C) IF
!                              ABS(Z) .GT. 1.
!
!                      C - +/-(A/R)
!
!                      S - +/-(B/R)
!
! MODULES REFERENCED BY GIVENS - NONE
!
! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT
!
!***********************************************************
!
      real aa , bb , r , u , v
!
! LOCAL PARAMETERS -
!
! AA,BB = LOCAL COPIES OF A AND B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   VARIABLES USED TO SCALE A AND B FOR COMPUTING R
!
      aa = a
      bb = b
      if ( abs(aa)>abs(bb) ) then
!
! ABS(A) .GT. ABS(B)
!
         u = aa + aa
         v = bb/u
         r = sqrt(.25+v*v)*u
         c = aa/r
         s = v*(c+c)
!
! NOTE THAT R HAS THE SIGN OF A, C .GT. 0, AND S HAS
!   SIGN(A)*SIGN(B)
!
         b = s
         a = r
         return
!
! ABS(A) .LE. ABS(B)
!
      elseif ( bb==0. ) then
!
! A = B = 0.
!
         c = 1.
         s = 0.
         goto 99999
      endif
      u = bb + bb
      v = aa/u
!
! STORE R IN A
!
      a = sqrt(.25+v*v)*u
      s = bb/a
      c = v*(s+s)
!
! NOTE THAT R HAS THE SIGN OF B, S .GT. 0, AND C HAS
!   SIGN(A)*SIGN(B)
!
      b = 1.
      if ( c/=0. ) b = 1./c
      return
99999 end subroutine givens
!*==GRADG.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine gradg(n,x,y,z,w,iadj,iend,eps,nit,grad,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - W
!A INPUT  - IADJ
!A INPUT  - IEND
!A INPUT  - EPS
!A OUTPUT - NIT
!A OUTPUT - GRAD
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       APLYRT   CONSTR
! called by   AA0001   AA0003
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  A11      A12      A22      ALFA     CX       CY       D        DG1      DG2      DGK      DGMAX    G1       G2       G3       INDF     INDL     INDX     ITER     K        MAXIT    NB       NN       R1       R2       SINAL    SX       SY
!             T        TOL      WK       XK       XNB      XS       YK       YNB      YS       ZK       ZNB
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iadj(*) , iend(n) , nit , ier
      real x(n) , y(n) , z(n) , w(n) , eps , grad(3,n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N NODES ON THE UNIT SPHERE WITH
! ASSOCIATED DATA VALUES, THIS ROUTINE USES A GLOBAL METHOD
! TO COMPUTE ESTIMATED GRADIENTS AT THE NODES.  THE METHOD
! CONSISTS OF MINIMIZING A QUADRATIC FUNCTIONAL Q(G) OVER
! THE N-VECTOR G OF GRADIENTS WHERE Q APPROXIMATES THE LIN-
! EARIZED CURVATURE OF AN INTERPOLANT F OVER THE TRIANGULA-
! TION.  THE RESTRICTION OF F TO AN ARC OF THE TRIANGULATION
! IS TAKEN TO BE THE HERMITE CUBIC (WITH RESPECT TO ARC-
! LENGTH) INTERPOLANT OF THE DATA VALUES AND TANGENTIAL
! GRADIENT COMPONENTS AT THE ENDPOINTS OF THE ARC, AND Q IS
! THE SUM OF THE LINEARIZED CURVATURES OF F ALONG THE ARCS
! -- THE INTEGRALS OVER THE ARCS OF D2F(A)**2 WHERE D2F(A)
! IS THE SECOND DERIVATIVE OF F WITH RESPECT TO ARC-LENGTH
! A.
!   SINCE THE GRADIENT AT NODE K LIES IN THE PLANE TANGENT
! TO THE SPHERE SURFACE AT K, IT IS DEFINED BY TWO COMPO-
! NENTS -- ITS X AND Y COMPONENTS IN THE COORDINATE SYSTEM
! OBTAINED BY ROTATING K TO THE NORTH POLE.  THUS, THE MIN-
! IMIZATION PROBLEM CORRESPONDS TO AN ORDER 2N SYMMETRIC
! POSITIVE-DEFINITE SPARSE LINEAR SYSTEM WHICH IS SOLVED BY
! THE BLOCK GAUSS-SEIDEL METHOD WITH 2 BY 2 BLOCKS.
!   AN ALTERNATIVE METHOD, SUBROUTINE GRADL, COMPUTES A
! LOCAL APPROXIMATION TO THE GRADIENT AT A SINGLE NODE AND,
! WHILE SLIGHTLY LESS EFFICIENT, WAS FOUND TO BE GENERALLY
! MORE ACCURATE WHEN THE NODAL DISTRIBUTION IS VERY DENSE,
! VARIES GREATLY, OR DOES NOT COVER THE SPHERE.  GRADG, ON
! THE OTHER HAND, WAS FOUND TO BE SLIGHTLY MORE ACCURATE ON
! A SOMEWHAT UNIFORM DISTRIBUTION OF 514 NODES.
!
! INPUT PARAMETERS - N - NUMBER OF NODES.  N .GE. 3.
!
!                X,Y,Z - CARTESIAN COORDINATES OF THE NODES.
!                        X(I)**2 + Y(I)**2 + Z(I)**2 = 1 FOR
!                        I = 1,...,N.
!
!                    W - DATA VALUES AT THE NODES.  W(I) IS
!                        ASSOCIATED WITH (X(I),Y(I),Z(I)).
!
!            IADJ,IEND - DATA STRUCTURE DEFINING THE TRIAN-
!                        GULATION.  SEE SUBROUTINE TRMESH.
!
!                  EPS - NONNEGATIVE CONVERGENCE CRITERION.
!                        THE METHOD IS TERMINATED WHEN THE
!                        MAXIMUM CHANGE IN A GRADIENT COMPO-
!                        NENT BETWEEN ITERATIONS IS AT MOST
!                        EPS.  EPS = 1.E-3 IS SUFFICIENT FOR
!                        EFFECTIVE CONVERGENCE.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                  NIT - MAXIMUM NUMBER OF GAUSS-SEIDEL
!                        ITERATIONS TO BE APPLIED.  THIS
!                        MAXIMUM WILL LIKELY BE ACHIEVED IF
!                        EPS IS SMALLER THAN THE MACHINE
!                        PRECISION.  OPTIMAL EFFICIENCY WAS
!                        ACHIEVED IN TESTING WITH EPS = 0
!                        AND NIT = 5 OR 6.
!
!                 GRAD - INITIAL SOLUTION ESTIMATES (ZERO
!                        VECTORS ARE SUFFICIENT).  GRAD(I,J)
!                        CONTAINS COMPONENT I OF THE GRADI-
!                        ENT AT NODE J FOR I = 1,2,3 (X,Y,Z)
!                        AND J = 1,...,N.  GRAD( ,J) MUST BE
!                        ORTHOGONAL TO NODE J -- GRAD(1,J)*
!                        X(J) + GRAD(2,J)*Y(J) + GRAD(3,J)*
!                        Z(J) = 0.
!
! OUTPUT PARAMETERS - NIT - NUMBER OF GAUSS-SEIDEL ITERA-
!                           TIONS EMPLOYED.
!
!                    GRAD - ESTIMATED GRADIENTS.  SEE THE
!                           DESCRIPTION UNDER INPUT PARAME-
!                           TERS.  GRAD IS NOT CHANGED IF
!                           IER = 2.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF THE CONVERGENCE CRI-
!                                   TERION WAS ACHIEVED.
!                           IER = 1 IF CONVERGENCE WAS NOT
!                                   ACHIEVED WITHIN NIT
!                                   ITERATIONS.
!                           IER = 2 IF N OR EPS IS OUT OF
!                                   RANGE OR NIT .LT. 0 ON
!                                   INPUT.
!
! MODULES REFERENCED BY GRADG - CONSTR, APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY GRADG - ATAN, SQRT, AMAX1,
!                                       ABS
!
!***********************************************************
!
      integer nn , maxit , iter , k , indf , indl , indx , nb
      real tol , dgmax , xk , yk , zk , wk , g1 , g2 , g3 , cx , sx ,   &
     &     cy , sy , a11 , a12 , a22 , r1 , r2 , xnb , ynb , znb ,      &
     &     alfa , xs , ys , sinal , d , t , dg1 , dg2 , dgk(3)
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! MAXIT =       INPUT VALUE OF NIT
! ITER =        NUMBER OF ITERATIONS USED
! K =           DO-LOOP AND NODE INDEX
! INDF,INDL =   IADJ INDICES OF THE FIRST AND LAST NEIGHBORS
!                 OF K
! INDX =        IADJ INDEX IN THE RANGE INDF,...,INDL
! NB =          NEIGHBOR OF K
! TOL =         LOCAL COPY OF EPS
! DGMAX =       MAXIMUM CHANGE IN A GRADIENT COMPONENT BE-
!                 TWEEN ITERATIONS
! XK,YK,ZK,WK = X(K), Y(K), Z(K), W(K)
! G1,G2,G3 =    COMPONENTS OF GRAD( ,K)
! CX,SX,CY,SY = COMPONENTS OF A ROTATION MAPPING K TO THE
!                 NORTH POLE (0,0,1)
! A11,A12,A22 = MATRIX COMPONENTS OF THE 2 BY 2 BLOCK A*DG
!                 = R WHERE A IS SYMMETRIC, (DG1,DG2,0) IS
!                 THE CHANGE IN THE GRADIENT AT K, AND R IS
!                 THE RESIDUAL
! R1,R2 =       COMPONENTS OF THE RESIDUAL -- DERIVATIVES OF
!                 Q WITH RESPECT TO THE COMPONENTS OF THE
!                 GRADIENT AT NODE K
! XNB,YNB,ZNB = COORDINATES OF NODE NB IN THE ROTATED COOR-
!                 DINATE SYSTEM
! ALFA =        ARC-LENGTH BETWEEN NODES K AND NB
! XS,YS =       XNB**2, YNB**2
! SINAL =       SIN(ALFA) -- MAGNITUDE OF THE VECTOR CROSS
!                 PRODUCT K X NB
! D =           ALFA*SINAL**2 -- FACTOR IN THE 2 BY 2 SYSTEM
! T =           TEMPORARY STORAGE AND FACTOR OF R1 AND R2
! DG1,DG2 =     SOLUTION OF THE 2 BY 2 SYSTEM -- FIRST 2
!                 COMPONENTS OF DGK IN THE ROTATED COORDI-
!                 NATE SYSTEM
! DGK =         CHANGE IN GRAD( ,K) FROM THE PREVIOUS ESTI-
!                 MATE
!
      nn = n
      tol = eps
      maxit = nit
!
! ERROR CHECKS AND INITIALIZATION
!
      if ( nn<3 .or. tol<0. .or. maxit<0 ) then
!
! PARAMETER OUT OF RANGE
!
         nit = 0
         ier = 2
         goto 99999
      else
         iter = 0
!
! TOP OF ITERATION LOOP
!
         do while ( iter/=maxit )
            dgmax = 0.
            indl = 0
            do k = 1 , nn
               xk = x(k)
               yk = y(k)
               zk = z(k)
               wk = w(k)
               g1 = grad(1,k)
               g2 = grad(2,k)
               g3 = grad(3,k)
!
!   CONSTRUCT THE ROTATION MAPPING NODE K TO THE NORTH POLE
!
               call constr(x(k),y(k),z(k),cx,sx,cy,sy)
!
!   INITIALIZE COMPONENTS OF THE 2 BY 2 SYSTEM
!
               a11 = 0.
               a12 = 0.
               a22 = 0.
               r1 = 0.
               r2 = 0.
!
!   LOOP ON NEIGHBORS NB OF K
!
               indf = indl + 1
               indl = iend(k)
               do indx = indf , indl
                  nb = iadj(indx)
                  if ( nb/=0 ) then
!
!   COMPUTE THE COORDINATES OF NB IN THE ROTATED SYSTEM
!
                     t = sx*y(nb) + cx*z(nb)
                     ynb = cx*y(nb) - sx*z(nb)
                     znb = sy*x(nb) + cy*t
                     xnb = cy*x(nb) - sy*t
!
!   COMPUTE ARC-LENGTH ALFA BETWEN NB AND K, SINAL =
!     SIN(ALFA), AND D = ALFA*SIN(ALFA)**2
!
                     alfa = 2.*atan(sqrt((1.-znb)/(1.+znb)))
                     xs = xnb*xnb
                     ys = ynb*ynb
                     sinal = sqrt(xs+ys)
                     d = alfa*(xs+ys)
!
!   UPDATE THE SYSTEM COMPONENTS FOR NODE NB
!
                     a11 = a11 + xs/d
                     a12 = a12 + xnb*ynb/d
                     a22 = a22 + ys/d
                     t = 1.5*(w(nb)-wk)/(alfa*alfa*sinal)               &
     &                   + ((grad(1,nb)*xk+grad(2,nb)*yk+grad(3,nb)*zk) &
     &                   /2.-(g1*x(nb)+g2*y(nb)+g3*z(nb)))/d
                     r1 = r1 + t*xnb
                     r2 = r2 + t*ynb
                  endif
               enddo
!
!   SOLVE THE 2 BY 2 SYSTEM AND UPDATE DGMAX
!
               dg2 = (a11*r2-a12*r1)/(a11*a22-a12*a12)
               dg1 = (r1-a12*dg2)/a11
               dgmax = amax1(dgmax,abs(dg1),abs(dg2))
!
!   ROTATE (DG1,DG2,0) BACK TO THE ORIGINAL COORDINATE
!     SYSTEM AND UPDATE GRAD( ,K)
!
               call aplyrt(dg1,dg2,cx,sx,cy,sy,dgk)
               grad(1,k) = g1 + dgk(1)
               grad(2,k) = g2 + dgk(2)
               grad(3,k) = g3 + dgk(3)
            enddo
!
!   INCREMENT ITER AND TEST FOR CONVERGENCE
!
            iter = iter + 1
            if ( dgmax<=tol ) then
!
! METHOD CONVERGED
!
               nit = iter
               ier = 0
               return
            endif
         enddo
      endif
!
! METHOD FAILED TO CONVERGE WITHIN NIT ITERATIONS
!
      ier = 1
      return
99999 end subroutine gradg
!*==GRADL.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine gradl(n,k,x,y,z,w,iadj,iend,g,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - K
!A PASSED - X
!A PASSED - Y
!A PASSED - Z
!A INPUT  - W
!A PASSED - IADJ
!A PASSED - IEND
!A PASSED - G
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       APLYR    APLYRT   CONSTR   GETNP    GIVENS   ROTATE   SETUP
! called by   AA0001   AA0003   INTRC1   INTRC1   UNIF     UNIF
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  A        AV       AVSQ     C        CX       CY       DF       DMIN     DTOL     DX       DY       I        IERR     IM1      IP1      J        JP1      KK       L        LM1      LMAX     LMIN     LMN      LMX      LNP      NN       NP
!             NPTS     RF       RIN      RTOL     S        SF       SUM      SX       SY       WK       WT       XP       YP       ZP
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , k , iadj(*) , iend(n) , ier
      real x(n) , y(n) , z(n) , w(n) , g(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE WITH THEIR ASSOCIATED DATA VALUES W, THIS ROUTINE
! ESTIMATES A GRADIENT VECTOR AT NODE K AS FOLLOWS -- THE
! COORDINATE SYSTEM IS ROTATED SO THAT K BECOMES THE NORTH
! POLE, NODE K AND A SET OF NEARBY NODES ARE PROJECTED
! ORTHOGONALLY ONTO THE X-Y PLANE (IN THE NEW COORDINATE
! SYSTEM), A QUADRATIC IS FITTED IN A WEIGHTED LEAST-SQUARES
! SENSE TO THE DATA VALUES AT THE PROJECTED NODES SUCH THAT
! THE VALUE (ASSOCIATED WITH K) AT (0,0) IS INTERPOLATED, X-
! AND Y-PARTIAL DERIVATIVE ESTIMATES DX AND DY ARE COMPUTED
! BY DIFFERENTIATING THE QUADRATIC AT (0,0), AND THE ESTI-
! MATED GRADIENT G IS OBTAINED BY ROTATING (DX,DY,0) BACK TO
! THE ORIGINAL COORDINATE SYSTEM.  NOTE THAT G LIES IN THE
! PLANE TANGENT TO THE SPHERE AT NODE K, I.E. G IS ORTHOGO-
! NAL TO THE UNIT VECTOR REPRESENTED BY NODE K.  A MARQUARDT
! STABILIZATION FACTOR IS USED IF NECESSARY TO ENSURE A
! WELL-CONDITIONED LEAST SQUARES SYSTEM, AND A UNIQUE SOLU-
! TION EXISTS UNLESS THE NODES ARE COLLINEAR.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES IN THE TRIANGU-
!                           LATION.  N .GE. 7.
!
!                       K - NODE AT WHICH THE GRADIENT IS
!                           SOUGHT.  1 .LE. K .LE. N.
!
!                   X,Y,Z - CARTESIAN COORDINATES OF THE
!                           NODES.
!
!                       W - DATA VALUES AT THE NODES.  W(I)
!                           IS ASSOCIATED WITH (X(I),Y(I),
!                           Z(I)) FOR I = 1,...,N.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           IN THE TRIANGULATION.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS FOR EACH NODE.
!
! IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -    G - X-, Y-, AND Z-COMPONENTS (IN
!                            THAT ORDER) OF THE ESTIMATED
!                            GRADIENT AT NODE K UNLESS
!                            IER .LT. 0.
!
!                      IER - ERROR INDICATOR
!                            IER .GE. 6 IF NO ERRORS WERE
!                                       ENCOUNTERED.  IER
!                                       CONTAINS THE NUMBER
!                                       OF NODES (INCLUDING
!                                       K) USED IN THE LEAST
!                                       SQUARES FIT.
!                            IER = -1 IF N OR K IS OUT OF
!                                     RANGE.
!                            IER = -2 IF THE LEAST SQUARES
!                                     SYSTEM HAS NO UNIQUE
!                                     SOLUTION DUE TO DUP-
!                                     LICATE OR COLLINEAR
!                                     NODES.
!
! MODULES REFERENCED BY GRADL - GETNP, CONSTR, APLYR,
!                               SETUP, GIVENS, ROTATE,
!                               APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY GRADL - MIN0, FLOAT, SQRT,
!                                       AMIN1, ABS
!
!***********************************************************
!
      integer nn , kk , lmn , lmx , lmin , lmax , lm1 , lnp , npts(30) ,&
     &        ierr , np , i , j , im1 , ip1 , jp1 , l
      real wk , sum , df , rf , rtol , avsq , av , rin , cx , sx , cy , &
     &     sy , xp , yp , zp , wt , a(6,6) , c , s , dmin , dtol , sf , &
     &     dx , dy
      data lmn/10/
      data lmx/30/ , rtol/1.E-6/ , dtol/.01/ , sf/1./
!
! LOCAL PARAMETERS -
!
! NN,KK =     LOCAL COPIES OF N AND K
! LMN,LMX =   MINIMUM AND MAXIMUM VALUES OF LNP FOR N
!               SUFFICIENTLY LARGE.  IN MOST CASES LMN-1
!               NODES ARE USED IN THE FIT.  7 .LE. LMN .LE.
!               LMX.
! LMIN,LMAX = MIN(LMN,N), MIN(LMX,N)
! LM1 =       LMIN-1
! LNP =       LENGTH OF NPTS OR LMAX+1
! NPTS =      ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!               NODES ORDERED BY ANGULAR DISTANCE FROM K.
!               NPTS(1)=K AND THE FIRST LNP-1 ELEMENTS OF
!               NPTS ARE USED IN THE LEAST SQUARES FIT.
!               UNLESS LNP = LMAX+1, NPTS(LNP) DETERMINES R
!               (SEE RIN).
! IERR =      ERROR FLAG FOR CALLS TO GETNP (NOT CHECKED)
! NP =        ELEMENT OF NPTS TO BE ADDED TO THE SYSTEM
! I,J =       LOOP INDICES
! IM1,IP1 =   I-1, I+1
! JP1 =       J+1
! L =         NUMBER OF COLUMNS OF A**T TO WHICH A ROTATION
!               IS APPLIED
! WK =        W(K) -- DATA VALUE AT NODE K
! SUM =       SUM OF SQUARED EUCLIDEAN DISTANCES (IN THE
!               ROTATED COORDINATE SYSTEM) BETWEEN THE
!               ORIGIN AND THE NODES USED IN THE LEAST
!               SQUARES FIT
! DF =        NEGATIVE Z-COMPONENT (IN THE ROTATED COORDI-
!               NATE SYSTEM) OF AN ELEMENT NP OF NPTS --
!               INCREASING FUNCTION OF THE ANGULAR DISTANCE
!               BETWEEN K AND NP.  DF LIES IN THE INTERVAL
!               (-1,1).
! RF =        VALUE OF DF ASSOCIATED WITH NPTS(LNP) UNLESS
!               LNP = LMAX+1 (SEE RIN)
! RTOL =      TOLERANCE FOR DETERMINING LNP (AND HENCE R) --
!               IF THE INCREASE IN DF BETWEEN TWO SUCCESSIVE
!               ELEMENTS OF NPTS IS LESS THAN RTOL, THEY ARE
!               TREATED AS BEING THE SAME DISTANCE FROM NODE
!               K AND AN ADDITIONAL NODE IS ADDED
! AVSQ =      AV*AV -- ACCUMULATED IN SUM
! AV =        ROOT-MEAN-SQUARE DISTANCE (IN THE ROTATED
!               COORDINATE SYSTEM) BETWEEN THE ORIGIN AND
!               THE NODES (OTHER THAN K) IN THE LEAST
!               SQUARES FIT.  THE FIRST 3 COLUMNS OF A**T
!               ARE SCALED BY 1/AVSQ, THE NEXT 2 BY 1/AV.
! RIN =       INVERSE OF A RADIUS OF INFLUENCE R WHICH
!               ENTERS INTO WT -- R = 1+RF UNLESS ALL ELE-
!               MENTS OF NPTS ARE USED IN THE FIT (LNP =
!               LMAX+1), IN WHICH CASE R IS THE DISTANCE
!               FUNCTION ASSOCIATED WITH SOME POINT MORE
!               DISTANT FROM K THAN NPTS(LMAX)
! CX,SX =     COMPONENTS OF A PLANE ROTATION ABOUT THE X-
!               AXIS WHICH, TOGETHER WITH CY AND SY, DEFINE
!               A MAPPING FROM NODE K TO THE NORTH POLE
!               (0,0,1)
! CY,SY =     COMPONENTS OF A PLANE ROTATION ABOUT THE Y-
!               AXIS
! XP,YP,ZP =  COORDINATES OF NP IN THE ROTATED COORDINATE
!               SYSTEM UNLESS ZP .LT. 0, IN WHICH CASE
!               (XP,YP,0) LIES ON THE EQUATOR
! WT =        WEIGHT FOR THE EQUATION CORRESPONDING TO NP --
!               WT = (R-D)/(R*D) = 1/D - RIN WHERE D = 1-ZP
!               IS ASSOCIATED WITH NP
! A =         TRANSPOSE OF THE (UPPER TRIANGLE OF THE) AUG-
!               MENTED REGRESSION MATRIX
! C,S =       COMPONENTS OF THE PLANE ROTATION USED TO
!               TRIANGULARIZE THE REGRESSION MATRIX
! DMIN =      MINIMUM OF THE MAGNITUDES OF THE DIAGONAL
!               ELEMENTS OF THE TRIANGULARIZED REGRESSION
!               MATRIX
! DTOL =      TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!               SYSTEM -- DMIN IS REQUIRED TO BE AT LEAST
!               DTOL
! SF =        MARQUARDT STABILIZATION FACTOR USED TO DAMP
!               OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!               PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!               IS ILL-CONDITIONED.  INCREASING SF RESULTS
!               IN MORE DAMPING (A MORE NEARLY LINEAR FIT).
! DX,DY =     X AND Y COMPONENTS OF THE ESTIMATED GRADIENT
!               IN THE ROTATED COORDINATE SYSTEM
!
      nn = n
      kk = k
      wk = w(kk)
!
! CHECK FOR ERRORS AND INITIALIZE LMIN, LMAX
!
      if ( nn<7 .or. kk<1 .or. kk>nn ) then
!
! N OR K IS OUT OF RANGE
!
         ier = -1
         return
      else
         lmin = min0(lmn,nn)
         lmax = min0(lmx,nn)
!
! COMPUTE NPTS, LNP, AVSQ, AV, AND R.
!   SET NPTS TO THE CLOSEST LMIN-1 NODES TO K.  DF CONTAINS
!   THE NEGATIVE Z-COMPONENT (IN THE ROTATED COORDINATE
!   SYSTEM) OF THE NEW NODE ON RETURN FROM GETNP.
!
         sum = 0.
         npts(1) = kk
         lm1 = lmin - 1
         do lnp = 2 , lm1
            call getnp(x,y,z,iadj,iend,lnp,npts,df,ierr)
            sum = sum + 1. - df*df
         enddo
!
!   ADD ADDITIONAL NODES TO NPTS UNTIL THE INCREASE IN
!     R = 1+RF IS AT LEAST RTOL.
!
         do lnp = lmin , lmax
            call getnp(x,y,z,iadj,iend,lnp,npts,rf,ierr)
            if ( rf-df>=rtol ) goto 50
            sum = sum + 1. - rf*rf
         enddo
!
!   USE ALL LMAX NODES IN THE LEAST SQUARES FIT.  R IS
!     ARBITRARILY INCREASED BY 5 PER CENT.
!
         rf = 1.05*rf + .05
         lnp = lmax + 1
!
!   THERE ARE LNP-2 EQUATIONS CORRESPONDING TO NODES
!     NPTS(2),...,NPTS(LNP-1).
!
 50      avsq = sum/float(lnp-2)
         av = sqrt(avsq)
         rin = 1./(1.+rf)
!
! CONSTRUCT THE ROTATION
!
         call constr(x(kk),y(kk),z(kk),cx,sx,cy,sy)
!
! SET UP THE FIRST 5 EQUATIONS OF THE AUGMENTED REGRESSION
!   MATRIX (TRANSPOSED) AS THE COLUMNS OF A, AND ZERO OUT
!   THE LOWER TRIANGLE (UPPER TRIANGLE OF A) WITH GIVENS
!   ROTATIONS
!
         do i = 1 , 5
            np = npts(i+1)
            call aplyr(x(np),y(np),z(np),cx,sx,cy,sy,xp,yp,zp)
            wt = 1./(1.-zp) - rin
            call setup(xp,yp,w(np),wk,av,avsq,wt,a(1,i))
            if ( i/=1 ) then
               im1 = i - 1
               do j = 1 , im1
                  jp1 = j + 1
                  l = 6 - j
                  call givens(a(j,j),a(j,i),c,s)
                  call rotate(l,c,s,a(jp1,j),a(jp1,i))
               enddo
            endif
         enddo
!
! ADD THE ADDITIONAL EQUATIONS TO THE SYSTEM USING
!   THE LAST COLUMN OF A -- I .LE. LNP.
!
         i = 7
         do
            if ( i==lnp ) then
!
! TEST THE SYSTEM FOR ILL-CONDITIONING
!
               dmin = amin1(abs(a(1,1)),abs(a(2,2)),abs(a(3,3)),        &
     &                abs(a(4,4)),abs(a(5,5)))
               if ( dmin>=dtol ) exit
               if ( lnp>lmax ) then
!
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS --ADD
!   MULTIPLES OF THE FIRST THREE UNIT VECTORS TO THE FIRST
!   THREE EQUATIONS.
!
                  do i = 1 , 3
                     a(i,6) = sf
                     ip1 = i + 1
                     do j = ip1 , 6
                        a(j,6) = 0.
                     enddo
                     do j = i , 5
                        jp1 = j + 1
                        l = 6 - j
                        call givens(a(j,j),a(j,6),c,s)
                        call rotate(l,c,s,a(jp1,j),a(jp1,6))
                     enddo
                  enddo
!
! TEST THE LINEAR PORTION OF THE STABILIZED SYSTEM FOR
!   ILL-CONDITIONING
!
                  dmin = amin1(abs(a(4,4)),abs(a(5,5)))
                  if ( dmin>=dtol ) exit
!
! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES
!
                  ier = -2
                  goto 99999
               else
!
! ADD ANOTHER NODE TO THE SYSTEM AND INCREASE R --
!   I .EQ. LNP
!
                  lnp = lnp + 1
                  if ( lnp<=lmax ) call getnp(x,y,z,iadj,iend,lnp,npts, &
     &                 rf,ierr)
                  rin = 1./(1.05*(1.+rf))
               endif
            else
               np = npts(i)
               call aplyr(x(np),y(np),z(np),cx,sx,cy,sy,xp,yp,zp)
               wt = 1./(1.-zp) - rin
               call setup(xp,yp,w(np),wk,av,avsq,wt,a(1,6))
               do j = 1 , 5
                  jp1 = j + 1
                  l = 6 - j
                  call givens(a(j,j),a(j,6),c,s)
                  call rotate(l,c,s,a(jp1,j),a(jp1,6))
               enddo
               i = i + 1
            endif
         enddo
      endif
!
! SOLVE THE 2 BY 2 TRIANGULAR SYSTEM FOR THE ESTIMATED
!   PARTIAL DERIVATIVES
!
      dy = a(6,5)/a(5,5)
      dx = (a(6,4)-a(5,4)*dy)/a(4,4)/av
      dy = dy/av
!
! ROTATE THE GRADIENT (DX,DY,0) BACK INTO THE ORIGINAL
!   COORDINATE SYSTEM
!
      call aplyrt(dx,dy,cx,sx,cy,sy,g)
      ier = lnp - 1
      return
99999 end subroutine gradl
!*==INDEX.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      integer function indexr(nvertx,nabor,iadj,iend)
!*--********************************************************************
!A INPUT  - NVERTX
!A INPUT  - NABOR
!A INPUT  - IADJ
!A INPUT  - IEND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   ADNODE   ADNODE   DELETE   DELETE   SWAP     SWAP
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  INDX     NB
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer nvertx , nabor , iadj(*) , iend(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION RETURNS THE INDEX OF NABOR IN THE
! ADJACENCY LIST FOR NVERTX.
!
! INPUT PARAMETERS - NVERTX - NODE WHOSE ADJACENCY LIST IS
!                             TO BE SEARCHED.
!
!                     NABOR - NODE WHOSE INDEX IS TO BE
!                             RETURNED.  NABOR MUST BE
!                             CONNECTED TO NVERTX.
!
!                      IADJ - SET OF ADJACENCY LISTS.
!
!                      IEND - POINTERS TO THE ENDS OF
!                             ADJACENCY LISTS IN IADJ.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER -  INDEX - IADJ(INDEX) = NABOR.
!
! MODULES REFERENCED BY INDEX - NONE
!
!***********************************************************
!
      integer nb , indx
!
! LOCAL PARAMETERS -
!
! NB =   LOCAL COPY OF NABOR
! INDX = INDEX FOR IADJ
!
      nb = nabor
!
! INITIALIZATION
!
      indx = iend(nvertx) + 1
      do
!
! SEARCH THE LIST OF NVERTX NEIGHBORS FOR NB
!
         indx = indx - 1
         if ( iadj(indx)==nb ) then
!
            indexr = indx
            exit
         endif
      enddo
      end function indexr
!*==INTADD.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine intadd(kk,i1,i2,i3,iadj,iend)
!*--********************************************************************
!A INPUT  - KK
!A INPUT  - I1
!A INPUT  - I2
!A INPUT  - I3
!A OUTPUT - IADJ
!A OUTPUT - IEND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       SHIFTD
! called by   ADNODE   ADNODE
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IMAX     IMIN     INDX     IP1      IP2      IP3      ITEMP    K        KM1      N        N1       N2       NF       NFT      NL
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer kk , i1 , i2 , i3 , iadj(*) , iend(kk)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS AN INTERIOR NODE TO A TRIANGULATION
! OF A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND
! ARE UPDATED WITH THE INSERTION OF NODE KK IN THE TRIANGLE
! WHOSE VERTICES ARE I1, I2, AND I3.
!
! INPUT PARAMETERS -        KK - INDEX OF NODE TO BE
!                                INSERTED.  KK .GE. 4.
!
!                     I1,I2,I3 - INDICES OF THE VERTICES OF
!                                A TRIANGLE CONTAINING NODE
!                                KK - IN COUNTERCLOCKWISE
!                                ORDER.
!
!                         IADJ - SET OF ADJACENCY LISTS
!                                OF NODES IN THE MESH.
!
!                         IEND - POINTERS TO THE ENDS OF
!                                ADJACENCY LISTS IN IADJ FOR
!                                EACH NODE IN THE MESH.
!
!   IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN
! THE VERTICES I1, I2, AND I3.  I1,I2,I3 MAY BE DETERMINED
! BY TRFIND.
!
! KK, I1, I2, AND I3 ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ADDITION
!                                 OF NODE KK AS THE LAST
!                                 ENTRY.  NODE KK WILL BE
!                                 CONNECTED TO NODES I1, I2,
!                                 AND I3.  NO OPTIMIZATION
!                                 OF THE MESH IS PERFORMED.
!
! MODULE REFERENCED BY INTADD - SHIFTD
!
! INTRINSIC FUNCTION CALLED BY INTADD - MOD
!
!***********************************************************
!
      integer k , km1 , n(3) , nft(3) , ip1 , ip2 , ip3 , indx , nf ,   &
     &        nl , n1 , n2 , imin , imax , i , itemp
!
! LOCAL PARAMETERS -
!
! K =           LOCAL COPY OF KK
! KM1 =         K - 1
! N =           VECTOR CONTAINING I1, I2, I3
! NFT =         POINTERS TO THE TOPS OF THE 3 SETS OF IADJ
!                 ELEMENTS TO BE SHIFTED DOWNWARD
! IP1,IP2,IP3 = PERMUTATION INDICES FOR N AND NFT
! INDX =        INDEX FOR IADJ AND N
! NF,NL =       INDICES OF FIRST AND LAST ENTRIES IN IADJ
!                 TO BE SHIFTED DOWN
! N1,N2 =       FIRST 2 VERTICES OF A NEW TRIANGLE --
!                 (N1,N2,KK)
! IMIN,IMAX =   BOUNDS ON DO-LOOP INDEX -- FIRST AND LAST
!                 ELEMENTS OF IEND TO BE INCREMENTED
! I =           DO-LOOP INDEX
! ITEMP =       TEMPORARY STORAGE LOCATION
!
      k = kk
!
! INITIALIZATION
!
      n(1) = i1
      n(2) = i2
      n(3) = i3
!
! SET UP NFT
!
      do i = 1 , 3
         n1 = n(i)
         indx = mod(i,3) + 1
         n2 = n(indx)
         indx = iend(n1) + 1
         do
!
! FIND THE INDEX OF N2 AS A NEIGHBOR OF N1
!
            indx = indx - 1
            if ( iadj(indx)==n2 ) then
               nft(i) = indx + 1
               exit
            endif
         enddo
      enddo
!
! ORDER THE VERTICES BY DECREASING MAGNITUDE -
!   N(IP(I+1)) PRECEDES N(IP(I)) IN IEND FOR
!   I = 1,2.
!
      ip1 = 1
      ip2 = 2
      ip3 = 3
      if ( n(2)>n(1) ) then
         ip1 = 2
         ip2 = 1
      endif
      if ( n(3)>n(ip1) ) then
         ip3 = ip1
         ip1 = 3
      endif
      if ( n(ip3)>n(ip2) ) then
         itemp = ip2
         ip2 = ip3
         ip3 = itemp
      endif
!
! ADD NODE K TO THE ADJACENCY LISTS OF EACH VERTEX AND
!   UPDATE IEND.  FOR EACH VERTEX, A SET OF IADJ ELEMENTS
!   IS SHIFTED DOWNWARD AND K IS INSERTED.  SHIFTING STARTS
!   AT THE END OF THE ARRAY.
!
      km1 = k - 1
      nl = iend(km1)
      nf = nft(ip1)
      if ( nf<=nl ) call shiftd(nf,nl,3,iadj)
      iadj(nf+2) = k
      imin = n(ip1)
      imax = km1
      do i = imin , imax
         iend(i) = iend(i) + 3
      enddo
!
      nl = nf - 1
      nf = nft(ip2)
      call shiftd(nf,nl,2,iadj)
      iadj(nf+1) = k
      imax = imin - 1
      imin = n(ip2)
      do i = imin , imax
         iend(i) = iend(i) + 2
      enddo
!
      nl = nf - 1
      nf = nft(ip3)
      call shiftd(nf,nl,1,iadj)
      iadj(nf) = k
      imax = imin - 1
      imin = n(ip3)
      do i = imin , imax
         iend(i) = iend(i) + 1
      enddo
!
! ADD NODE K TO IEND AND ITS NEIGHBORS TO IADJ
!
      indx = iend(km1)
      iend(k) = indx + 3
      do i = 1 , 3
         indx = indx + 1
         iadj(indx) = n(i)
      enddo
      end subroutine intadd
!*==INTRC0.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine intrc0(n,plat,plon,x,y,z,w,iadj,iend,ist,pw,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - PLAT
!A INPUT  - PLON
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - W
!A INPUT  - IADJ
!A INPUT  - IEND
!A OUTPUT - IST
!A OUTPUT - PW
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       TRFIND
! called by   AA0001   AA0003
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  B1       B2       B3       I1       I2       I3       INDX     N1       N2       P        P1       P2       P3       PTN1     PTN2     S12      SUM
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iadj(*) , iend(n) , ist , ier
      real plat , plon , x(n) , y(n) , z(n) , w(n) , pw
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE COMPUTES THE VALUE AT A POINT P OF A CONTINUOUS
! FUNCTION WHICH INTERPOLATES THE DATA VALUES.  THE INTERP-
! OLATORY FUNCTION IS LINEAR ON EACH UNDERLYING TRIANGLE
! (PLANAR TRIANGLE WITH THE SAME VERTICES AS A SPHERICAL
! TRIANGLE).  IF P IS NOT CONTAINED IN A TRIANGLE, AN EX-
! TRAPOLATED VALUE IS TAKEN TO BE THE INTERPOLATED VALUE AT
! THE NEAREST POINT OF THE TRIANGULATION BOUNDARY.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE TRIANGU-
!                            LATION.  N .GE. 3.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!                IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                            CREATED BY SUBROUTINE TRMESH.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
!                      PW - VALUE OF THE INTERPOLATORY
!                           FUNCTION AT P IF IER .LE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = 1 IF EXTRAPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -1 IF N .LT. 3 OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!                           IER = -3 IF P IS NOT IN A TRI-
!                                    ANGLE AND THE ANGLE BE-
!                                    TWEEN P AND THE NEAREST
!                                    BOUNDARY POINT IS AT
!                                    LEAST 90 DEGREES.
!
! MODULES REFERENCED BY INTRC0 - TRFIND
!
! INTRINSIC FUNCTIONS CALLED BY INTRC0 - COS, SIN
!
!***********************************************************
!
      integer i1 , i2 , i3 , n1 , n2 , indx
      real p(3) , p1(3) , p2(3) , p3(3) , b1 , b2 , b3 , sum , s12 ,    &
     &     ptn1 , ptn2
!
! LOCAL PARAMETERS -
!
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! N1,N2 =    ENDPOINTS OF A BOUNDARY ARC WHICH IS VISIBLE
!              FROM P WHEN P IS NOT CONTAINED IN A TRIANGLE
! INDX =     IADJ INDEX OF N1 AS A NEIGHBOR OF N2
! P =        CARTESIAN COORDINATES OF P
! P1,P2,P3 = CARTESIAN COORDINATES OF I1, I2, AND I3
! B1,B2,B3 = BARYCENTRIC COORDINATES OF THE CENTRAL PROJEC-
!              TION OF P ONTO THE UNDERLYING PLANAR TRIANGLE
!              OR (B1 AND B2) PROJECTION OF Q ONTO THE
!              UNDERLYING LINE SEGMENT N1-N2 WHEN P IS
!              EXTERIOR -- UNNORMALIZED COORDINATES ARE
!              COMPUTED BY TRFIND WHEN P IS IN A TRIANGLE
! SUM =      QUANTITY USED TO NORMALIZE THE BARYCENTRIC
!              COORDINATES
! S12 =      SCALAR PRODUCT (N1,N2)
! PTN1 =     SCALAR PRODUCT (P,N1)
! PTN2 =     SCALAR PRODUCT (P,N2)
!
      if ( n<3 .or. ist<1 .or. ist>n ) then
!
! N OR IST OUT OF RANGE
!
         ier = -1
         return
      else
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
         p(1) = cos(plat)*cos(plon)
         p(2) = cos(plat)*sin(plon)
         p(3) = sin(plat)
!
! FIND THE VERTEX INDICES OF A TRIANGLE CONTAINING P
!
         call trfind(ist,p,x,y,z,iadj,iend,b1,b2,b3,i1,i2,i3)
         if ( i1==0 ) then
!
! COLLINEAR NODES
!
            ier = -2
            return
         else
            ist = i1
            if ( i3==0 ) then
!
! P IS EXTERIOR TO THE TRIANGULATION, AND I1 AND I2 ARE
!   BOUNDARY NODES WHICH ARE VISIBLE FROM P.  SET PW TO THE
!   INTERPOLATED VALUE AT Q WHERE Q IS THE CLOSEST BOUNDARY
!   POINT TO P.
!
! TRAVERSE THE BOUNDARY STARTING FROM THE RIGHTMOST VISIBLE
!   NODE I1.
!
               n1 = i1
               ptn1 = p(1)*x(n1) + p(2)*y(n1) + p(3)*z(n1)
               if ( i1==i2 ) then
                  do
!
! ALL BOUNDARY NODES ARE VISIBLE FROM P.  FIND A BOUNDARY
!   ARC N1->N2 SUCH THAT P LEFT (N2 X N1)->N1
!
! COUNTERCLOCKWISE BOUNDARY TRAVERSAL --
!   SET N2 TO THE FIRST NEIGHBOR OF N1
!
                     indx = 1
                     if ( n1/=1 ) indx = iend(n1-1) + 1
                     n2 = iadj(indx)
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N2), AND COMPUTE
!   B2 = DET(P,N1,N2 X N1)
!
                     s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
                     ptn2 = p(1)*x(n2) + p(2)*y(n2) + p(3)*z(n2)
                     b2 = ptn2 - s12*ptn1
                     if ( b2<=0. ) exit
!
! P RIGHT (N2 X N1)->N1 -- ITERATE
!
                     n1 = n2
                     i1 = n1
                     ptn1 = ptn2
                  enddo
               endif
               do
!
! P LEFT (N2 X N1)->N1 WHERE N2 IS THE FIRST NEIGHBOR OF P1
!   CLOCKWISE BOUNDARY TRAVERSAL --
!
                  n2 = n1
                  ptn2 = ptn1
!
! SET N1 TO THE LAST NEIGHBOR OF N2 AND TEST FOR TERMINATION
!
                  indx = iend(n2) - 1
                  n1 = iadj(indx)
                  if ( n1==i1 ) then
!
! THE ANGULAR DISTANCE BETWEEN P AND THE CLOSEST BOUNDARY
!   POINT TO P IS AT LEAST 90 DEGREES
!
                     ier = -3
                     goto 99999
                  else
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N1)
!
                     s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
                     ptn1 = p(1)*x(n1) + p(2)*y(n1) + p(3)*z(n1)
!
! COMPUTE B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q)
!
                     b2 = ptn2 - s12*ptn1
                     if ( b2>0. ) then
!
! COMPUTE B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q)
!
                        b1 = ptn1 - s12*ptn2
                        if ( b1>0. ) exit
!
! Q = N2
!
                        pw = w(n2)
                        ier = 1
                        return
                     endif
                  endif
               enddo
            else
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE
!   VERTEX COORDINATES IN LOCAL VARIABLES
!
               p1(1) = x(i1)
               p1(2) = y(i1)
               p1(3) = z(i1)
               p2(1) = x(i2)
               p2(2) = y(i2)
               p2(3) = z(i2)
               p3(1) = x(i3)
               p3(2) = y(i3)
               p3(3) = z(i3)
!
! NORMALIZE THE BARYCENTRIC COORDINATES
!
               sum = b1 + b2 + b3
               b1 = b1/sum
               b2 = b2/sum
               b3 = b3/sum
               pw = b1*w(i1) + b2*w(i2) + b3*w(i3)
               ier = 0
               return
            endif
         endif
      endif
!
! P STRICTLY LEFT (N2 X N1)->N2 AND P STRICTLY LEFT
!   N1->(N2 X N1).  THUS Q LIES ON THE INTERIOR OF N1->N2.
!   NORMALIZE THE COORDINATES AND COMPUTE PW.
!
      sum = b1 + b2
      pw = (b1*w(n1)+b2*w(n2))/sum
      ier = 1
      return
99999 end subroutine intrc0
!*==INTRC1.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine intrc1(n,plat,plon,x,y,z,w,iadj,iend,iflag,grad,ist,pw,&
     &                  ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - PLAT
!A INPUT  - PLON
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - W
!A INPUT  - IADJ
!A INPUT  - IEND
!A INPUT  - IFLAG
!A INPUT  - GRAD
!A OUTPUT - IST
!A OUTPUT - PW
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ARCINT   ARCLEN   GRADL    TRFIND   WVAL
! called by   UNIF     UNIF
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  A        B1       B2       B3       DUM      G1       G2       G3       GQ       GQN      I        I1       I2       I3       IERR     INDX     N1       N2       NN       P        P1       P2       P3       PTGQ     PTN1     PTN2     Q
!             QNORM    S12      SUM      W1       W2       W3       WQ
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      !real arclen
!*** End of declarations inserted by SPAG
      integer n , iadj(*) , iend(n) , iflag , ist , ier
      real plat , plon , x(n) , y(n) , z(n) , w(n) , grad(3,n) , pw
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF NODES ON THE UNIT
! SPHERE, ALONG WITH DATA VALUES AT THE NODES, THIS SUB-
! ROUTINE CONSTRUCTS THE VALUE AT A POINT P OF A ONCE CON-
! TINUOUSLY DIFFERENTIABLE FUNCTION WHICH INTERPOLATES THE
! DATA VALUES.  IF THE TRIANGULATION DOES NOT COVER THE
! ENTIRE SPHERE, THE SURFACE IS EXTENDED CONTINUOUSLY BEYOND
! THE BOUNDARY ALLOWING EXTRAPOLATION.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES.  N .GE. 3 AND
!                            N .GE. 7 IF IFLAG = 0.
!
!                PLAT,PLON - LATITUDE AND LONGITUDE OF P IN
!                            RADIANS.
!
!                    X,Y,Z - VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF THE NODES.
!
!                        W - VECTOR CONTAINING DATA VALUES
!                            AT THE NODES.  W(I) IS ASSOCI-
!                            ATED WITH (X(I),Y(I),Z(I)) FOR
!                            I = 1,...,N.
!
!                IADJ,IEND - DATA STRUCTURE REPRESENTING THE
!                            TRIANGULATION.  SEE SUBROUTINE
!                            TRMESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF INTRC1 IS TO PRO-
!                                      VIDE ESTIMATED GRAD-
!                                      IENTS (FROM GRADL).
!                                      N .GE. 7 IN THIS
!                                      CASE.
!                            IFLAG = 1 IF GRADIENTS ARE PRO-
!                                      VIDED IN GRAD.  THIS
!                                      IS MORE EFFICIENT IF
!                                      INTRC1 IS TO BE
!                                      CALLED SEVERAL TIMES.
!
!                     GRAD - ARRAY DIMENSIONED 3 BY N WHOSE
!                            I-TH COLUMN CONTAINS AN ESTI-
!                            MATED GRADIENT AT NODE I IF
!                            IFLAG = 1 (SEE SUBROUTINE
!                            GRADL).  GRAD MAY BE A DUMMY
!                            VARIABLE (NOT USED) IF IFLAG
!                            = 0.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING P.  1 .LE. IST .LE. N.
!                            THE OUTPUT VALUE OF IST FROM A
!                            PREVIOUS CALL MAY BE A GOOD
!                            CHOICE.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING P (OR
!                           NEAREST P) UNLESS IER = -1 OR
!                           IER = -2.
!
! OUTPUT PARAMETERS -  PW - INTERPOLATED VALUE AT P IF
!                           IER .GE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF INTERPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = 1 IF EXTRAPOLATION WAS
!                                   PERFORMED SUCCESSFULLY.
!                           IER = -1 IF N, IFLAG, OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!                           IER = -3 IF THE ANGULAR DISTANCE
!                                    BETWEEN P AND THE NEAR-
!                                    EST POINT OF THE TRIAN-
!                                    GULATION IS AT LEAST 90
!                                    DEGREES.
!
! MODULES REFERENCED BY INTRC1 - TRFIND, WVAL, ARCINT,
!                                ARCLEN,
!             (AND OPTIONALLY)   GRADL, GETNP, CONSTR,
!                                APLYR, SETUP, GIVENS,
!                                ROTATE, APLYRT
!
! INTRINSIC FUNCTIONS CALLED BY INTRC1 - COS, SIN, SQRT
!
!***********************************************************
!
      integer nn , i1 , i2 , i3 , i , ierr , n1 , n2 , indx
      real p(3) , p1(3) , p2(3) , p3(3) , w1 , w2 , w3 , g1(3) , g2(3) ,&
     &     g3(3) , b1 , b2 , b3 , sum , dum(3) , s12 , ptn1 , ptn2 ,    &
     &     q(3) , qnorm , wq , gq(3) , a , ptgq , gqn
!
! LOCAL PARAMETERS -
!
! NN =       LOCAL COPY OF N
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! I =        DO-LOOP INDEX
! IERR =     ERROR FLAG FOR CALLS TO GRADL
! N1,N2 =    INDICES OF THE ENDPOINTS OF A BOUNDARY ARC WHEN
!              P IS EXTERIOR (NOT CONTAINED IN A TRIANGLE)
! INDX =     IADJ INDEX OF N2 AS A NEIGHBOR OF N1 OR VICE-
!              VERSA
! P =        CARTESIAN COORDINATES OF P
! P1,P2,P3 = CARTESIAN COORDINATES OF THE VERTICES I1, I2,
!              AND I3, OR (P1 AND P2) COORDINATES OF N1 AND
!              N2 IF P IS EXTERIOR
! W1,W2,W3 = DATA VALUES ASSOCIATED WITH I1, I2, AND I3, OR
!              (W1 AND W2) VALUES ASSOCIATED WITH N1 AND
!              N2 IF P IS EXTERIOR
! G1,G2,G3 = GRADIENTS AT I1, I2, AND I3, OR (G1 AND G2) AT
!              N1 AND N2
! B1,B2,B3 = BARYCENTRIC COORDINATES OF THE CENTRAL PROJEC-
!              TION OF P ONTO THE UNDERLYING PLANAR TRIANGLE
!              OR (B1 AND B2) PROJECTION OF Q ONTO THE
!              UNDERLYING LINE SEGMENT N1-N2 WHEN P IS
!              EXTERIOR -- UNNORMALIZED COORDINATES ARE
!              COMPUTED BY TRFIND WHEN P IS IN A TRIANGLE
! SUM =      QUANTITY USED TO NORMALIZE THE BARYCENTRIC
!              COORDINATES
! DUM =      DUMMY PARAMETER FOR WVAL AND ARCINT
! S12 =      SCALAR PRODUCT (N1,N2) -- FACTOR OF B1 AND B2
! PTN1 =     SCALAR PRODUCT (P,N1) -- FACTOR OF B1 AND B2
! PTN2 =     SCALAR PRODUCT (P,N2) -- FACTOR OF B1 AND B2
! Q =        CLOSEST BOUNDARY POINT TO P WHEN P IS EXTERIOR
! QNORM =    FACTOR USED TO NORMALIZE Q
! WQ,GQ =    INTERPOLATED VALUE AND GRADIENT AT Q
! A =        ANGULAR SEPARATION BETWEEN P AND Q
! PTGQ =     SCALAR PRODUCT (P,GQ) -- FACTOR OF THE COMPONENT
!              OF GQ IN THE DIRECTION Q->P
! GQN =      NEGATIVE OF THE COMPONENT OF GQ IN THE DIRECTION
!              Q->P
!
      nn = n
      if ( nn<3 .or. (iflag/=1 .and. nn<7) .or. iflag<0 .or.            &
     &     iflag>1 .or. ist<1 .or. ist>nn ) then
!
! N, IFLAG, OR IST OUT OF RANGE
!
         ier = -1
         return
      else
!
! TRANSFORM (PLAT,PLON) TO CARTESIAN COORDINATES
!
         p(1) = cos(plat)*cos(plon)
         p(2) = cos(plat)*sin(plon)
         p(3) = sin(plat)
!
! LOCATE P WITH RESPECT TO THE TRIANGULATION
!
         call trfind(ist,p,x,y,z,iadj,iend,b1,b2,b3,i1,i2,i3)
         if ( i1/=0 ) then
            ist = i1
            if ( i3==0 ) then
!
! P IS EXTERIOR TO THE TRIANGULATION, AND I1 AND I2 ARE
!   BOUNDARY NODES WHICH ARE VISIBLE FROM P.  EXTRAPOLATE TO
!   P BY LINEAR (WITH RESPECT TO ARC-LENGTH) INTERPOLATION
!   OF THE VALUE AND DIRECTIONAL DERIVATIVE (GRADIENT COMP-
!   ONENT IN THE DIRECTION Q->P) OF THE INTERPOLATORY
!   SURFACE AT Q WHERE Q IS THE CLOSEST BOUNDARY POINT TO P.
!
! DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM I1
!
               n1 = i1
               ptn1 = p(1)*x(n1) + p(2)*y(n1) + p(3)*z(n1)
               if ( i1==i2 ) then
                  do
!
! ALL BOUNDARY NODES ARE VISIBLE FROM P.  FIND A BOUNDARY
!   ARC N1->N2 SUCH THAT P LEFT (N2 X N1)->N1
!
! COUNTERCLOCKWISE BOUNDARY TRAVERSAL --
!   SET N2 TO THE FIRST NEIGHBOR OF N1
!
                     indx = 1
                     if ( n1/=1 ) indx = iend(n1-1) + 1
                     n2 = iadj(indx)
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N2), AND COMPUTE
!   B2 = DET(P,N1,N2 X N1)
!
                     s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
                     ptn2 = p(1)*x(n2) + p(2)*y(n2) + p(3)*z(n2)
                     b2 = ptn2 - s12*ptn1
                     if ( b2<=0. ) exit
!
! P RIGHT (N2 X N1)->N1 -- ITERATE
!
                     n1 = n2
                     i1 = n1
                     ptn1 = ptn2
                  enddo
               endif
               do
!
! P LEFT (N2 X N1)->N1 WHERE N2 IS THE FIRST NEIGHBOR OF N1
!   CLOCKWISE BOUNDARY TRAVERSAL --
!
                  n2 = n1
                  ptn2 = ptn1
!
! SET N1 TO THE LAST NEIGHBOR OF N2 AND TEST FOR TERMINATION
!
                  indx = iend(n2) - 1
                  n1 = iadj(indx)
                  if ( n1==i1 ) then
!
! THE DISTANCE BETWEEN P AND THE CLOSEST BOUNDARY POINT
!   IS AT LEAST 90 DEGREES
!
                     ier = -3
                     goto 99999
                  else
!
! COMPUTE INNER PRODUCTS (N1,N2) AND (P,N1)
!
                     s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
                     ptn1 = p(1)*x(n1) + p(2)*y(n1) + p(3)*z(n1)
!
! COMPUTE B2 = DET(P,N1,N2 X N1) = DET(Q,N1,N2 X N1)*(P,Q)
!
                     b2 = ptn2 - s12*ptn1
                     if ( b2>0. ) then
!
! COMPUTE B1 = DET(P,N2 X N1,N2) = DET(Q,N2 X N1,N2)*(P,Q)
!
                        b1 = ptn1 - s12*ptn2
                        if ( b1>0. ) then
!
! P STRICTLY LEFT (N2 X N1)->N2 AND P STRICTLY LEFT
!   N1->(N2 X N1).  THUS Q LIES ON THE INTERIOR OF N1->N2.
!   STORE DATA VALUES AND COORDINATES OF N1 AND N2 IN
!   LOCAL VARIABLES
!
                           w1 = w(n1)
                           w2 = w(n2)
                           p1(1) = x(n1)
                           p1(2) = y(n1)
                           p1(3) = z(n1)
                           p2(1) = x(n2)
                           p2(2) = y(n2)
                           p2(3) = z(n2)
!
! COMPUTE THE CENTRAL PROJECTION OF Q ONTO P2-P1 AND
!   NORMALIZE TO OBTAIN Q
!
                           qnorm = 0.
                           do i = 1 , 3
                              q(i) = b1*p1(i) + b2*p2(i)
                              qnorm = qnorm + q(i)*q(i)
                           enddo
                           qnorm = sqrt(qnorm)
                           do i = 1 , 3
                              q(i) = q(i)/qnorm
                           enddo
!
! STORE OR COMPUTE GRADIENTS AT N1 AND N2
!
                           if ( iflag/=1 ) then
                              call gradl(nn,n1,x,y,z,w,iadj,iend,g1,    &
     &                           ierr)
                              if ( ierr<0 ) exit
                              call gradl(nn,n2,x,y,z,w,iadj,iend,g2,    &
     &                           ierr)
                              if ( ierr<0 ) exit
                           else
                              do i = 1 , 3
                                 g1(i) = grad(i,n1)
                                 g2(i) = grad(i,n2)
                              enddo
                           endif
!
! COMPUTE AN INTERPOLATED VALUE AND NORMAL GRADIENT-
!   COMPONENT AT Q
!
                           call arcint(q,p1,p2,w1,w2,g1,g2,wq,dum,gqn)
!
! EXTRAPOLATE TO P -- THE NORMAL GRADIENT COMPONENT GQN IS
!   THE NEGATIVE OF THE COMPONENT IN THE DIRECTION Q->P.
!
                           pw = wq - gqn*arclen(q,p)
                           ier = 1
                           return
                        else
!
! Q = N2.  STORE VALUE, COORDINATES, AND AND GRADIENT AT Q.
!
                           wq = w(n2)
                           q(1) = x(n2)
                           q(2) = y(n2)
                           q(3) = z(n2)
                           if ( iflag/=1 ) then
                              call gradl(nn,n2,x,y,z,w,iadj,iend,gq,    &
     &                           ierr)
                              if ( ierr<0 ) exit
                           else
                              do i = 1 , 3
                                 gq(i) = grad(i,n2)
                              enddo
                           endif
!
! EXTRAPOLATE TO P -- PW = WQ + A*(GQ,Q X (PXQ)/SIN(A))
!   WHERE A IS THE ANGULAR SEPARATION BETWEEN Q AND P,
!   AND SIN(A) IS THE MAGNITUDE OF P X Q.
!
                           a = arclen(q,p)
                           ptgq = p(1)*gq(1) + p(2)*gq(2) + p(3)*gq(3)
                           pw = wq
                           if ( a/=0. ) pw = pw + ptgq*a/sin(a)
                           ier = 1
                           return
                        endif
                     endif
                  endif
               enddo
            else
!
! P IS CONTAINED IN THE TRIANGLE (I1,I2,I3).  STORE THE DATA
!   VALUES AND VERTEX COORDINATES IN LOCAL VARIABLES.
!
               w1 = w(i1)
               w2 = w(i2)
               w3 = w(i3)
               p1(1) = x(i1)
               p1(2) = y(i1)
               p1(3) = z(i1)
               p2(1) = x(i2)
               p2(2) = y(i2)
               p2(3) = z(i2)
               p3(1) = x(i3)
               p3(2) = y(i3)
               p3(3) = z(i3)
               if ( iflag/=1 ) then
!
! COMPUTE GRADIENT ESTIMATES AT THE VERTICES
!
                  call gradl(nn,i1,x,y,z,w,iadj,iend,g1,ierr)
                  if ( ierr<0 ) goto 100
                  call gradl(nn,i2,x,y,z,w,iadj,iend,g2,ierr)
                  if ( ierr<0 ) goto 100
                  call gradl(nn,i3,x,y,z,w,iadj,iend,g3,ierr)
                  if ( ierr<0 ) goto 100
               else
!
! GRADIENTS ARE USER-PROVIDED
!
                  do i = 1 , 3
                     g1(i) = grad(i,i1)
                     g2(i) = grad(i,i2)
                     g3(i) = grad(i,i3)
                  enddo
               endif
!
! NORMALIZE THE COORDINATES
!
               sum = b1 + b2 + b3
               b1 = b1/sum
               b2 = b2/sum
               b3 = b3/sum
               call wval(b1,b2,b3,p1,p2,p3,w1,w2,w3,g1,g2,g3,0,pw,dum)
               ier = 0
               return
            endif
         endif
      endif
!
! COLLINEAR NODES
!
 100  ier = -2
      return
99999 end subroutine intrc1
!*==PERMUT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine permut(n,ip,a)
!*--********************************************************************
!A INPUT  - N
!A OUTPUT - IP
!A OUTPUT - A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   REORDR   REORDR
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  IPJ      J        K        NN       TEMP
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , ip(n)
      real a(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE APPLIES A SET OF PERMUTATIONS TO A VECTOR.
!
! INPUT PARAMETERS -  N - LENGTH OF A AND IP.
!
!                    IP - VECTOR CONTAINING THE SEQUENCE OF
!                         INTEGERS 1,...,N PERMUTED IN THE
!                         SAME FASHION THAT A IS TO BE PER-
!                         MUTED.
!
!                     A - VECTOR TO BE PERMUTED.
!
! N AND IP ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -  A - REORDERED VECTOR REFLECTING THE
!                         PERMUTATIONS DEFINED BY IP.
!
! MODULES REFERENCED BY PERMUT - NONE
!
!***********************************************************
!
      integer nn , k , j , ipj
      real temp
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! K =    INDEX FOR IP AND FOR THE FIRST ELEMENT OF A IN A
!          PERMUTATION
! J =    INDEX FOR IP AND A, J .GE. K
! IPJ =  IP(J)
! TEMP = TEMPORARY STORAGE FOR A(K)
!
      nn = n
      if ( nn<2 ) return
      k = 1
!
! LOOP ON PERMUTATIONS
!
 100  j = k
      temp = a(k)
      do
!
! APPLY PERMUTATION TO A.  IP(J) IS MARKED (MADE NEGATIVE)
!   AS BEING INCLUDED IN THE PERMUTATION.
!
         ipj = ip(j)
         ip(j) = -ipj
         if ( ipj==k ) then
            a(j) = temp
            do
!
! SEARCH FOR AN UNMARKED ELEMENT OF IP
!
               k = k + 1
               if ( k>nn ) then
!
! ALL PERMUTATIONS HAVE BEEN APPLIED.  UNMARK IP.
!
                  do k = 1 , nn
                     ip(k) = -ip(k)
                  enddo
                  goto 99999
               elseif ( ip(k)>0 ) then
                  goto 100
               endif
            enddo
         else
            a(j) = a(ipj)
            j = ipj
         endif
      enddo
99999 end subroutine permut
!*==PRJCT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine prjct(px,py,pz,ox,oy,oz,ex,ey,ez,vx,vy,vz,isw,x,y,ier)
!*--********************************************************************
!A INPUT  - PX
!A INPUT  - PY
!A INPUT  - PZ
!A INPUT  - OX
!A INPUT  - OY
!A INPUT  - OZ
!A INPUT  - EX
!A INPUT  - EY
!A INPUT  - EZ
!A INPUT  - VX
!A INPUT  - VY
!A INPUT  - VZ
!A OUTPUT - ISW
!A OUTPUT - X
!A OUTPUT - Y
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   TRPLOT
! modifies    /PROCOM/ OES      XE       XH       XOE      XV       YE       YH       YOE      YV       ZE       ZH       ZOE      ZV
! uses value  /PROCOM/ OES      XE       XH       XOE      XV       YE       YH       YOE      YV       ZE       ZH       ZOE      ZV
! local vars  OEX      OEY      OEZ      S        SC       XEP      XW       YEP      YW       ZEP      ZW
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real oex , oey , oez
!*** End of declarations inserted by SPAG
      integer isw , ier
      real px , py , pz , ox , oy , oz , ex , ey , ez , vx , vy , vz ,  &
     &     x , y
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE PROJECTS A POINT P ONTO THE PLANE CON-
! TAINING POINT O WHOSE NORMAL IS THE LINE CONTAINING O AND
! E WHERE P, O, AND E ARE POINTS IN EUCLIDEAN 3-SPACE.
!
! INPUT PARAMETERS - PX,PY,PZ - CARTESIAN COORDINATES OF THE
!                               POINT P TO BE MAPPED ONTO
!                               THE PLANE.
!
!                    OX,OY,OZ - COORDINATES OF O (THE ORIGIN
!                               OF A COORDINATE SYSTEM IN
!                               THE PLANE OF PROJECTION).
!
!                    EX,EY,EZ - COORDINATES OF THE EYE-POSI-
!                               TION E DEFINING THE NORMAL
!                               TO THE PLANE AND THE LINE OF
!                               SIGHT FOR THE PROJECTION.  E
!                               MUST NOT COINCIDE WITH O OR
!                               P, AND THE ANGLE BETWEEN THE
!                               VECTORS O-E AND P-E MUST BE
!                               LESS THAN 90 DEGREES.  NOTE
!                               THAT E AND P MAY LIE ON OP-
!                               POSITE SIDES OF THE PLANE.
!
!                    VX,VY,VZ - COORDINATES OF A POINT V
!                               WHICH DEFINES THE POSITIVE
!                               Y-AXIS OF AN X-Y COORDINATE
!                               SYSTEM IN THE PLANE AS THE
!                               HALF-LINE CONTAINING O AND
!                               THE PROJECTION OF O+V ONTO
!                               THE PLANE.  THE POSITIVE X-
!                               AXIS HAS DIRECTION DEFINED
!                               BY THE CROSS PRODUCT V X
!                               E-O.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                         ISW - SWITCH WHICH MUST BE SET TO
!                               1 ON THE FIRST CALL AND WHEN
!                               THE VALUES OF OX,OY,OZ,EX,
!                               EY,EZ,VX,VY, OR VZ HAVE BEEN
!                               ALTERED SINCE A PREVIOUS
!                               CALL.  IF ISW .NE. 1, IT IS
!                               ASSUMED THAT ONLY THE COORD-
!                               INATES OF P HAVE CHANGED
!                               SINCE A PREVIOUS CALL.  PRE-
!                               VIOUSLY STORED QUANTITIES
!                               ARE USED FOR INCREASED EFF-
!                               ICIENCY IN THIS CASE.
!
! OUTPUT PARAMETERS - ISW - RESET TO 0 IF PRJCT WAS CALLED
!                           WITH ISW = 1, UNCHANGED OTHER-
!                           WISE.
!
!                     X,Y - COORDINATES OF THE POINT IN THE
!                           PLANE WHICH LIES ON THE LINE
!                           CONTAINING E AND P.  THE COORD-
!                           INATE SYSTEM IS DEFINED BY VX,
!                           VY, AND VZ.  X AND Y ARE NOT
!                           CHANGED IF IER .NE. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE
!                                   ENCOUNTERED.
!                           IER = 1 IF THE INNER PRODUCT OF
!                                   O-E WITH P-E IS NOT POS-
!                                   ITIVE IMPLYING THAT E IS
!                                   TOO CLOSE TO THE PLANE.
!                           IER = 2 IF O, E, AND O+V ARE
!                                   COLLINEAR.  SEE THE DES-
!                                   CRIPTION OF VX,VY,VZ.
!
! MODULES REFERENCED BY PRJCT - NONE
!
! INTRINSIC FUNCTION CALLED BY PRJCT - SQRT
!
!***********************************************************
!
      real xe , ye , ze , xoe , yoe , zoe , oes , s , xv , yv , zv ,    &
     &     sc , xh , yh , zh , xep , yep , zep , xw , yw , zw
      common /procom/ xe , ye , ze , xoe , yoe , zoe , oes , xv , yv ,  &
     &                zv , xh , yh , zh
!
! LOCAL PARAMETERS -
!
! XE,YE,ZE =    LOCAL COPIES OF EX, EY, EZ
! XOE,YOE,ZOE = COMPONENTS OF THE VECTOR OE FROM O TO E
! OES =         NORM SQUARED OF OE -- (OE,OE)
! S =           SCALE FACTOR FOR COMPUTING PROJECTIONS
! XV,YV,ZV =    COMPONENTS OF A UNIT VECTOR VN DEFINING THE
!                 POSITIVE Y-AXIS IN THE PLANE
! SC =          SCALE FACTOR FOR NORMALIZING VN AND HN
! XH,YH,ZH =    COMPONENTS OF A UNIT VECTOR HN DEFINING THE
!                 POSITIVE X-AXIS IN THE PLANE
! XEP,YEP,ZEP = COMPONENTS OF THE VECTOR EP FROM E TO P
! XW,YW,ZW =    COMPONENTS OF THE VECTOR W FROM O TO THE
!                 PROJECTION OF P ONTO THE PLANE
!
!   THE COMMON BLOCK IS USED TO SAVE LOCAL PARAMETERS
! BETWEEN CALLS.
!
      if ( isw==1 ) then
         isw = 0
!
! SET THE COORDINATES OF E TO LOCAL VARIABLES, COMPUTE
!   OE = E - O AND OES
!
         xe = ex
         ye = ey
         ze = ez
         xoe = xe - ox
         yoe = ye - oy
         zoe = ze - oz
         oes = xoe*xoe + yoe*yoe + zoe*zoe
         if ( oes==0. ) goto 100
!
! COMPUTE S = (OE,V)/OES AND VN = V - S*OE
!
         s = (xoe*vx+yoe*vy+zoe*vz)/oes
         xv = vx - s*oex
         yv = vy - s*oey
         zv = vz - s*oez
!
! NORMALIZE VN TO A UNIT VECTOR
!
         sc = xv*xv + yv*yv + zv*zv
         if ( sc==0. ) then
!
! O, E, AND O+V ARE COLLINEAR
!
            ier = 2
            goto 99999
         else
            sc = 1./sqrt(sc)
            xv = sc*xv
            yv = sc*yv
            zv = sc*zv
!
! COMPUTE HN = VN X OE (NORMALIZED)
!
            xh = yv*zoe - yoe*zv
            yh = xoe*zv - xv*zoe
            zh = xv*yoe - xoe*yv
            sc = 1./sqrt(xh*xh+yh*yh+zh*zh)
            xh = sc*xh
            yh = sc*yh
            zh = sc*zh
         endif
      endif
!
! COMPUTE EP = P - E, S = (OE,EP)/OES, AND W = OE - EP/S
!
      xep = px - xe
      yep = py - ye
      zep = pz - ze
      s = (xoe*xep+yoe*yep+zoe*zep)/oes
      if ( s<0. ) then
         xw = xoe - xep/s
         yw = yoe - yep/s
         zw = zoe - zep/s
!
! MAP W INTO X = (W,HN), Y = (W,VN)
!
         x = xw*xh + yw*yh + zw*zh
         y = xw*xv + yw*yv + zw*zv
         return
      endif
!
! (OE,EP) .GE. 0.
!
 100  ier = 1
      return
99999 end subroutine prjct
!*==QSORT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine qsort(n,x,ind)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A OUTPUT - IND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   REORDR   REORDR
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IJ       IL       INDX     IT       ITT      IU       J        K        L        M        R        T
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , ind(n)
      real x(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO
! SORT THE REAL ARRAY X INTO INCREASING ORDER.  THE ALGOR-
! ITHM IS AS FOLLOWS.  IND IS INITIALIZED TO THE ORDERED
! SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES ARE
! APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING
! A CENTRAL ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COM-
! PARED WITH T, AND INTERCHANGES ARE APPLIED AS NECESSARY SO
! THAT THE THREE VALUES ARE IN ASCENDING ORDER.  INTER-
! CHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS GREATER THAN
! T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER
! INDICES OF ONE OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS,
! AND THE PROCESS IS REPEATED ITERATIVELY ON THE OTHER
! PORTION.  WHEN A PORTION IS COMPLETELY SORTED, THE PROCESS
! BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.
!
! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.
!
!                      X - VECTOR OF LENGTH N TO BE SORTED.
!
!                    IND - VECTOR OF LENGTH .GE. N.
!
! N AND X ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N
!                          PERMUTED IN THE SAME FASHION AS X
!                          WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).
!
! MODULES REFERENCED BY QSORT - NONE
!
! INTRINSIC FUNCTIONS CALLED BY QSORT - IFIX, FLOAT
!
!***********************************************************
!
! NOTE -- IU AND IL MUST BE DIMENSIONED .GE. LOG(N) WHERE
!         LOG HAS BASE 2.
!
!***********************************************************
!
      integer iu(21) , il(21)
      integer m , i , j , k , l , ij , it , itt , indx
      real r , t
!
! LOCAL PARAMETERS -
!
! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X
!
      if ( n<=0 ) return
!
! INITIALIZE IND, M, I, J, AND R
!
      do i = 1 , n
         ind(i) = i
      enddo
      m = 1
      i = 1
      j = n
      r = .375
!
! TOP OF LOOP
!
 100  if ( i>=j ) goto 300
      if ( r>.5898437 ) then
         r = r - .21875
      else
         r = r + .0390625
      endif
!
! INITIALIZE K
!
 200  k = i
!
! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T
!
      ij = i + ifix(r*float(j-i))
      it = ind(ij)
      t = x(it)
!
! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T
!
      indx = ind(i)
      if ( x(indx)>t ) then
         ind(ij) = indx
         ind(i) = it
         it = indx
         t = x(it)
      endif
!
! INITIALIZE L
!
      l = j
!
! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T
!
      indx = ind(j)
      if ( x(indx)<t ) then
         ind(ij) = indx
         ind(j) = it
         it = indx
         t = x(it)
!
! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T
!
         indx = ind(i)
         if ( x(indx)>t ) then
            ind(ij) = indx
            ind(i) = it
            it = indx
            t = x(it)
         endif
      endif
      do
!
! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T
!
         l = l - 1
         indx = ind(l)
         if ( x(indx)<=t ) then
            do
!
! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS
!   NOT SMALLER THAN T
!
               k = k + 1
               indx = ind(k)
               if ( x(indx)>=t ) then
!
! IF K .LE. L, INTERCHANGE ELEMENTS K AND L
!
                  if ( k<=l ) then
!
! INTERCHANGE ELEMENTS K AND L
!
                     itt = ind(l)
                     ind(l) = ind(k)
                     ind(k) = itt
                     exit
                  else
!
! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED
!
                     if ( l-i<=j-k ) then
!
                        il(m) = k
                        iu(m) = j
                        j = l
                        m = m + 1
                     else
                        il(m) = i
                        iu(m) = l
                        i = k
                        m = m + 1
                     endif
                     goto 400
                  endif
               endif
            enddo
         endif
      enddo
!
! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY
!
 300  m = m - 1
      if ( m==0 ) return
      i = il(m)
      j = iu(m)
!
 400  if ( j-i>=11 ) goto 200
      if ( i==1 ) goto 100
      i = i - 1
      do
!
! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 .LE. I .LT. J AND
!   J-I .LT. 11.
!
         i = i + 1
         if ( i==j ) goto 300
         indx = ind(i+1)
         t = x(indx)
         it = indx
         indx = ind(i)
         if ( x(indx)>t ) then
            k = i
            do
!
               ind(k+1) = ind(k)
               k = k - 1
               indx = ind(k)
               if ( t>=x(indx) ) then
                  ind(k+1) = it
                  exit
               endif
            enddo
         endif
      enddo
      end subroutine qsort
!*==REORDR.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine reordr(n,iflag,a,b,c,d,ind)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - IFLAG
!A PASSED - A
!A PASSED - B
!A PASSED - C
!A PASSED - D
!A PASSED - IND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       PERMUT   QSORT
! called by   AA0001   AA0003   AA0004
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  NN       NV
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iflag , ind(n)
      real a(n) , b(n) , c(n) , d(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE USES A QUICK SORT TO REORDER THE REAL
! ARRAY A INTO INCREASING ORDER.  A RECORD OF THE PERMUTA-
! TIONS APPLIED TO A IS STORED IN IND, AND THESE PERMUTA-
! TIONS MAY BE APPLIED TO ONE, TWO, OR THREE ADDITIONAL
! VECTORS BY THIS ROUTINE.  ANY OTHER VECTOR V MAY BE PER-
! MUTED IN THE SAME FASHION BY CALLING SUBROUTINE PERMUT
! WITH N, IND, AND V AS PARAMETERS.
!   A SET OF NODES (X(I),Y(I),Z(I)) AND DATA VALUES W(I)
! MAY BE PREPROCESSED BY REORDR FOR INCREASED EFFICIENCY IN
! THE TRIANGULATION ROUTINE TRMESH.  THIS PREPROCESSING IS
! ALSO USEFUL FOR DETECTING DUPLICATE NODES.  NOTE THAT THE
! FIRST THREE NODES MUST NOT BE COLLINEAR ON INPUT TO
! TRMESH.  EITHER X, Y, OR Z MAY BE USED AS THE SORT KEY
! (ASSOCIATED WITH A).  IF THE NODAL COORDINATES ARE IN
! TERMS OF LATITUDE AND LONGITUDE, IT IS USUALLY ADVAN-
! TAGEOUS TO SORT ON LONGITUDE WITH THE SAME INTERCHANGES
! APPLIED TO LATITUDE AND THE DATA VALUES.  D WOULD BE A
! DUMMY PARAMETER IN THIS CASE.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES.
!
!                    IFLAG - NUMBER OF VECTORS TO BE PER-
!                            MUTED.
!                            IFLAG .LE. 0 IF A, B, C, AND
!                              D ARE TO REMAIN UNALTERED.
!                            IFLAG = 1 IF ONLY A IS TO BE
!                              PERMUTED.
!                            IFLAG = 2 IF A AND B ARE TO BE
!                              PERMUTED.
!                            IFLAG = 3 IF A, B, AND C ARE
!                              TO BE PERMUTED.
!                            IFLAG .GE. 4 IF A, B, C, AND D
!                              ARE TO BE PERMUTED.
!
!                  A,B,C,D - VECTORS OF LENGTH N TO BE
!                            SORTED (ON THE COMPONENTS OF A)
!                            OR DUMMY PARAMETERS, DEPENDING
!                            ON IFLAG.
!
!                      IND - VECTOR OF LENGTH .GE. N.
!
! N, IFLAG, AND ANY DUMMY PARAMETERS ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - A,B,C,D - SORTED OR UNALTERED VECTORS,
!                               DEPENDING ON IFLAG.
!
!                         IND - SEQUENCE OF INDICES 1,...,N
!                               PERMUTED IN THE SAME FASHION
!                               AS THE REAL VECTORS.  THUS
!                               THE ORDERING MAY BE APPLIED
!                               TO A VECTOR V AND STORED IN
!                               W BY SETTING W(I) =
!                               V(IND(I)), OR V MAY BE OVER-
!                               WRITTEN WITH THE ORDERING BY
!                               A CALL TO PERMUT.
!
! MODULES REFERENCED BY REORDR - QSORT, PERMUT
!
!***********************************************************
!
      integer nn , nv
!
! LOCAL PARAMETERS -
!
! NN = LOCAL COPY OF N
! NV = LOCAL COPY OF IFLAG
!
      nn = n
      nv = iflag
      call qsort(nn,a,ind)
      if ( nv<=0 ) return
      call permut(nn,ind,a)
      if ( nv==1 ) return
      call permut(nn,ind,b)
      if ( nv==2 ) return
      call permut(nn,ind,c)
      if ( nv==3 ) return
      call permut(nn,ind,d)
      end subroutine reordr
!*==ROTATE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine rotate(n,c,s,x,y)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - C
!A INPUT  - S
!A OUTPUT - X
!A OUTPUT - Y
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADL    GRADL    QSHEP2   QSHEP3
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        XI       YI
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n
      real c , s , x(n) , y(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!                                            ( C  S)
!   THIS ROUTINE APPLIES THE GIVENS ROTATION (     ) TO THE
!                                            (-S  C)
!               (X(1) ... X(N))
! 2 BY N MATRIX (             ).  THIS ROUTINE WAS TAKEN
!               (Y(1) ... Y(N))
! FROM LINPACK.
!
! INPUT PARAMETERS -   N - NUMBER OF COLUMNS TO BE ROTATED.
!
!                    C,S - ELEMENTS OF THE GIVENS ROTATION.
!                          THESE MAY BE DETERMINED BY
!                          SUBROUTINE GIVENS.
!
!                    X,Y - VECTORS OF LENGTH .GE. N
!                          CONTAINING THE 2-VECTORS TO BE
!                          ROTATED.
!
!   THE PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS
! ROUTINE.
!
! OUTPUT PARAMETERS - X,Y - ROTATED VECTORS
!
! MODULES REFERENCED BY ROTATE - NONE
!
!***********************************************************
!
      integer i
      real xi , yi
!
! LOCAL PARAMETERS -
!
! I =     DO-LOOP INDEX
! XI,YI = X(I), Y(I)
!
      if ( n<=0 .or. (c==1. .and. s==0.) ) return
      do i = 1 , n
         xi = x(i)
         yi = y(i)
         x(i) = c*xi + s*yi
         y(i) = -s*xi + c*yi
      enddo
      end subroutine rotate
!*==SETUP.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine setup(xi,yi,wi,wk,s1,s2,wt,row)
!*--********************************************************************
!A INPUT  - XI
!A INPUT  - YI
!A INPUT  - WI
!A INPUT  - WK
!A INPUT  - S1
!A INPUT  - S2
!A INPUT  - WT
!A OUTPUT - ROW
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   GRADL    GRADL
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  W1       W2
! uses PARAMs *** NONE ****
!*++********************************************************************
      real xi , yi , wi , wk , s1 , s2 , wt , row(6)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED
! REGRESSION MATRIX FOR A WEIGHTED LEAST SQUARES FIT OF A
! QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES WI
! WHERE Q(0,0) = WK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS)
! ARE SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LIN-
! EAR TERMS) ARE SCALED BY 1/S1.
!
! INPUT PARAMETERS - XI,YI - COORDINATES OF NODE I.
!
!                       WI - DATA VALUE AT NODE I.
!
!                       WK - DATA VALUE INTERPOLATED BY Q AT
!                            THE ORIGIN.
!
!                    S1,S2 - INVERSE SCALE FACTORS.
!
!                       WT - WEIGHT FACTOR CORRESPONDING TO
!                            THE I-TH EQUATION.
!
!                      ROW - VECTOR OF LENGTH 6.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - ROW - VECTOR CONTAINING A ROW OF THE
!                          AUGMENTED REGRESSION MATRIX.
!
! MODULES REFERENCED BY SETUP - NONE
!
!***********************************************************
!
      real w1 , w2
!
! LOCAL PARAMETERS -
!
! W1 = WEIGHTED SCALE FACTOR FOR THE LINEAR TERMS
! W2 = WEIGHTED SCALE FACTOR FOR THE QUADRATIC TERMS
!
      w1 = wt/s1
      w2 = wt/s2
      row(1) = xi*xi*w2
      row(2) = xi*yi*w2
      row(3) = yi*yi*w2
      row(4) = xi*w1
      row(5) = yi*w1
      row(6) = (wi-wk)*wt
      end subroutine setup
!*==SHIFTD.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine shiftd(nfrst,nlast,kk,iarr)
!*--********************************************************************
!A INPUT  - NFRST
!A INPUT  - NLAST
!A INPUT  - KK
!A OUTPUT - IARR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   BDYADD   BDYADD   DELETE   DELETE   INTADD   INTADD   SWAP     SWAP     TRMESH
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IBAK     IMAX     INC      INDX     K        NF       NL       NLP1     NS       NSL
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer nfrst , nlast , kk , iarr(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE SHIFTS A SET OF CONTIGUOUS ELEMENTS OF AN
! INTEGER ARRAY KK POSITIONS DOWNWARD (UPWARD IF KK .LT. 0).
! THE LOOPS ARE UNROLLED IN ORDER TO INCREASE EFFICIENCY.
!
! INPUT PARAMETERS - NFRST,NLAST - BOUNDS ON THE PORTION OF
!                                  IARR TO BE SHIFTED.  ALL
!                                  ELEMENTS BETWEEN AND
!                                  INCLUDING THE BOUNDS ARE
!                                  SHIFTED UNLESS NFRST .GT.
!                                  NLAST, IN WHICH CASE NO
!                                  SHIFT OCCURS.
!
!                             KK - NUMBER OF POSITIONS EACH
!                                  ELEMENT IS TO BE SHIFTED.
!                                  IF KK .LT. 0 SHIFT UP.
!                                  IF KK .GT. 0 SHIFT DOWN.
!
!                           IARR - INTEGER ARRAY OF LENGTH
!                                  .GE. NLAST + MAX(KK,0).
!
! NFRST, NLAST, AND KK ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -        IARR - SHIFTED ARRAY.
!
! MODULES REFERENCED BY SHIFTD - NONE
!
!***********************************************************
!
      integer inc , k , nf , nl , nlp1 , ns , nsl , i , ibak , indx ,   &
     &        imax
      data inc/5/
!
! LOCAL PARAMETERS -
!
! INC =  DO-LOOP INCREMENT (UNROLLING FACTOR) -- IF INC IS
!          CHANGED, STATEMENTS MUST BE ADDED TO OR DELETED
!          FROM THE DO-LOOPS
! K =    LOCAL COPY OF KK
! NF =   LOCAL COPY OF NFRST
! NL =   LOCAL COPY OF NLAST
! NLP1 = NL + 1
! NS =   NUMBER OF SHIFTS
! NSL =  NUMBER OF SHIFTS DONE IN UNROLLED DO-LOOP (MULTIPLE
!          OF INC)
! I =    DO-LOOP INDEX AND INDEX FOR IARR
! IBAK = INDEX FOR DOWNWARD SHIFT OF IARR
! INDX = INDEX FOR IARR
! IMAX = BOUND ON DO-LOOP INDEX
!
      k = kk
      nf = nfrst
      nl = nlast
      if ( nf>nl .or. k==0 ) return
      nlp1 = nl + 1
      ns = nlp1 - nf
      nsl = inc*(ns/inc)
      if ( k<0 ) then
!
! SHIFT UPWARD STARTING FROM THE TOP
!
         if ( nsl>0 ) then
            imax = nlp1 - inc
            do i = nf , imax , inc
               indx = i + k
               iarr(indx) = iarr(i)
               iarr(indx+1) = iarr(i+1)
               iarr(indx+2) = iarr(i+2)
               iarr(indx+3) = iarr(i+3)
               iarr(indx+4) = iarr(i+4)
            enddo
         endif
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
         i = nsl + nf
         do
            if ( i>nl ) return
            indx = i + k
            iarr(indx) = iarr(i)
            i = i + 1
         enddo
      else
!
! SHIFT DOWNWARD STARTING FROM THE BOTTOM
!
         if ( nsl>0 ) then
            do i = 1 , nsl , inc
               ibak = nlp1 - i
               indx = ibak + k
               iarr(indx) = iarr(ibak)
               iarr(indx-1) = iarr(ibak-1)
               iarr(indx-2) = iarr(ibak-2)
               iarr(indx-3) = iarr(ibak-3)
               iarr(indx-4) = iarr(ibak-4)
            enddo
         endif
!
! PERFORM THE REMAINING NS-NSL SHIFTS ONE AT A TIME
!
         ibak = nlp1 - nsl
         do
            if ( ibak<=nf ) return
            ibak = ibak - 1
            indx = ibak + k
            iarr(indx) = iarr(ibak)
         enddo
      endif
      end subroutine shiftd
!*==SWAP.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine swap(nin1,nin2,nout1,nout2,iadj,iend)
!*--********************************************************************
!A INPUT  - NIN1
!A INPUT  - NIN2
!A INPUT  - NOUT1
!A INPUT  - NOUT2
!A OUTPUT - IADJ
!A OUTPUT - IEND
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       INDEX    SHIFTD
! called by   ADNODE   ADNODE   EDGE     EDGE
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IMAX     IMIN     IN       IO       IP1      IP2      J        K        NF       NL
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      !integer index
!*** End of declarations inserted by SPAG
      integer nin1 , nin2 , nout1 , nout2 , iadj(*) , iend(*)
      !external index
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE SWAPS THE DIAGONALS IN A CONVEX QUADRI-
! LATERAL.
!
! INPUT PARAMETERS -  NIN1,NIN2,NOUT1,NOUT2 - NODAL INDICES
!                            OF A PAIR OF ADJACENT TRIANGLES
!                            WHICH FORM A CONVEX QUADRILAT-
!                            ERAL.  NOUT1 AND NOUT2 ARE CON-
!                            NECTED BY AN ARC WHICH IS TO BE
!                            REPLACED BY THE ARC NIN1-NIN2.
!                            (NIN1,NOUT1,NOUT2) MUST BE TRI-
!                            ANGLE VERTICES IN COUNTERCLOCK-
!                            WISE ORDER.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                            (SEE SUBROUTINE TRMESH).
!
! OUTPUT PARAMETERS - IADJ,IEND - UPDATED WITH THE ARC
!                                 REPLACEMENT.
!
! MODULES REFERENCED BY SWAP - INDEX, SHIFTD
!
!***********************************************************
!
      integer in(2) , io(2) , ip1 , ip2 , j , k , nf , nl , i , imin ,  &
     &        imax
!
! LOCAL PARAMETERS -
!
! IN =        NIN1 AND NIN2 ORDERED BY INCREASING MAGNITUDE
!               (THE NEIGHBORS OF IN(1) PRECEDE THOSE OF
!               IN(2) IN IADJ)
! IO =        NOUT1 AND NOUT2 IN INCREASING ORDER
! IP1,IP2 =   PERMUTATION OF (1,2) SUCH THAT IO(IP1)
!               PRECEDES IO(IP2) AS A NEIGHBOR OF IN(1)
! J,K =       PERMUTATION OF (1,2) USED AS INDICES OF IN
!               AND IO
! NF,NL =     IADJ INDICES BOUNDARY A PORTION OF THE ARRAY
!               TO BE SHIFTED
! I =         IEND INDEX
! IMIN,IMAX = BOUNDS ON THE PORTION OF IEND TO BE INCRE-
!               MENTED OR DECREMENTED
!
      in(1) = nin1
      in(2) = nin2
      io(1) = nout1
      io(2) = nout2
      ip1 = 1
!
! ORDER THE INDICES SO THAT IN(1) .LT. IN(2) AND IO(1) .LT.
!   IO(2), AND CHOOSE IP1 AND IP2 SUCH THAT (IN(1),IO(IP1),
!   IO(IP2)) FORMS A TRIANGLE.
!
      if ( in(1)>=in(2) ) then
         in(1) = in(2)
         in(2) = nin1
         ip1 = 2
      endif
      if ( io(1)>=io(2) ) then
         io(1) = io(2)
         io(2) = nout1
         ip1 = 3 - ip1
      endif
      ip2 = 3 - ip1
      if ( io(2)<in(1) ) then
!
! THE VERTICES ARE ORDERED (IO(1),IO(2),IN(1),IN(2)).
!   DELETE IO(2) BY SHIFTING UP BY 1
!
         nf = 1 + indexr(io(1),io(2),iadj,iend)
         nl = -1 + indexr(io(2),io(1),iadj,iend)
         if ( nf<=nl ) call shiftd(nf,nl,-1,iadj)
         imin = io(1)
         imax = io(2) - 1
         do i = imin , imax
            iend(i) = iend(i) - 1
         enddo
!
!   DELETE IO(1) BY SHIFTING UP BY 2 AND INSERT IN(2)
!
         nf = nl + 2
         nl = -1 + indexr(in(1),io(ip2),iadj,iend)
         if ( nf<=nl ) call shiftd(nf,nl,-2,iadj)
         iadj(nl-1) = in(2)
         imin = io(2)
         imax = in(1) - 1
         do i = imin , imax
            iend(i) = iend(i) - 2
         enddo
!
!   SHIFT UP BY 1 AND INSERT IN(1)
!
         nf = nl + 1
         nl = -1 + indexr(in(2),io(ip1),iadj,iend)
         call shiftd(nf,nl,-1,iadj)
         iadj(nl) = in(1)
         imin = in(1)
         imax = in(2) - 1
         do i = imin , imax
            iend(i) = iend(i) - 1
         enddo
         return
      elseif ( in(2)<io(1) ) then
!
! THE VERTICES ARE ORDERED (IN(1),IN(2),IO(1),IO(2)).
!   DELETE IO(1) BY SHIFTING DOWN BY 1
!
         nf = 1 + indexr(io(1),io(2),iadj,iend)
         nl = -1 + indexr(io(2),io(1),iadj,iend)
         if ( nf<=nl ) call shiftd(nf,nl,1,iadj)
         imin = io(1)
         imax = io(2) - 1
         do i = imin , imax
            iend(i) = iend(i) + 1
         enddo
!
!   DELETE IO(2) BY SHIFTING DOWN BY 2 AND INSERT IN(1)
!
         nl = nf - 2
         nf = 1 + indexr(in(2),io(ip2),iadj,iend)
         if ( nf<=nl ) call shiftd(nf,nl,2,iadj)
         iadj(nf+1) = in(1)
         imin = in(2)
         imax = io(1) - 1
         do i = imin , imax
            iend(i) = iend(i) + 2
         enddo
!
!   SHIFT DOWN BY 1 AND INSERT IN(2)
!
         nl = nf - 1
         nf = 1 + indexr(in(1),io(ip1),iadj,iend)
         call shiftd(nf,nl,1,iadj)
         iadj(nf) = in(2)
         imin = in(1)
         imax = in(2) - 1
         do i = imin , imax
            iend(i) = iend(i) + 1
         enddo
         goto 99999
      endif
!
! IN(1) AND IO(1) PRECEDE IN(2) AND IO(2).  FOR (J,K) =
!   (1,2) AND (2,1), DELETE IO(K) AS A NEIGHBOR OF IO(J)
!   BY SHIFTING A PORTION OF IADJ EITHER UP OR DOWN AND
!   AND INSERT IN(K) AS A NEIGHBOR OF IN(J).
!
      do j = 1 , 2
         k = 3 - j
         if ( in(j)>io(j) ) then
!
!   THE NEIGHBORS OF IO(J) PRECEDE THOSE OF IN(J) -- SHIFT
!     UP BY 1
!
            nf = 1 + indexr(io(j),io(k),iadj,iend)
            nl = -1 + indexr(in(j),io(ip2),iadj,iend)
            if ( nf<=nl ) call shiftd(nf,nl,-1,iadj)
            iadj(nl) = in(k)
            imin = io(j)
            imax = in(j) - 1
            do i = imin , imax
               iend(i) = iend(i) - 1
            enddo
         else
!
!   THE NEIGHBORS OF IN(J) PRECEDE THOSE OF IO(J) -- SHIFT
!     DOWN BY 1
!
            nf = 1 + indexr(in(j),io(ip1),iadj,iend)
            nl = -1 + indexr(io(j),io(k),iadj,iend)
            if ( nf<=nl ) call shiftd(nf,nl,1,iadj)
            iadj(nf) = in(k)
            imin = in(j)
            imax = io(j) - 1
            do i = imin , imax
               iend(i) = iend(i) + 1
            enddo
         endif
!
!   REVERSE (IP1,IP2) FOR (J,K) = (2,1)
!
         ip1 = ip2
         ip2 = 3 - ip1
      enddo
      return
99999 end subroutine swap
!*==SWPTST.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      logical function swptst(n1,n2,n3,n4,x,y,z)
!*--********************************************************************
!A INPUT  - N1
!A INPUT  - N2
!A INPUT  - N3
!A INPUT  - N4
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   ADNODE   ADNODE   ARCTST   ARCTST   EDGE     EDGE
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DX1      DX2      DX3      DY1      DY2      DY3      DZ1      DZ2      DZ3      X4       Y4       Z4
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n1 , n2 , n3 , n4
      real x(*) , y(*) , z(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION DECIDES WHETHER OR NOT TO REPLACE A
! DIAGONAL ARC IN A QUADRILATERAL WITH THE OTHER DIAGONAL.
! THE DECISION WILL BE POSITIVE (SWPTST = .TRUE.) IF AND
! ONLY IF N4 LIES ABOVE THE PLANE (IN THE HALF-SPACE NOT
! CONTAINING THE ORIGIN) DEFINED BY (N1,N2,N3), OR EQUIV-
! ALENTLY, IF THE PROJECTION OF N4 ONTO THIS PLANE IS
! INTERIOR TO THE CIRCUMCIRCLE OF (N1,N2,N3).  THE DECISION
! WILL BE NEGATIVE IF THE QUADRILATERAL IS NOT STRICTLY
! CONVEX.
!
! INPUT PARAMETERS - N1,N2,N3,N4 - INDICES OF THE FOUR NODES
!                      DEFINING THE QUADRILATERAL.  N1 AND
!                      N2 ARE CURRENTLY CONNECTED BY A DIAG-
!                      ONAL ARC WHICH SHOULD BE REPLACED BY
!                      AN ARC CONNECTING N3 TO N4 IF THE
!                      DECISION IS MADE TO SWAP.  N1, N2, N3
!                      MUST BE IN COUNTERCLOCKWISE ORDER.
!
!              X,Y,Z - VECTORS OF NODAL COORDINATES.  (X(I),
!                      Y(I),Z(I)) ARE THE COORDINATES OF
!                      NODE I FOR I = N1, N2, N3, AND N4.
!
! NONE OF THE INPUT PARAMETERS ARE ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -  SWPTST - .TRUE. IFF THE ARC CONNECTING
!                              N1 AND N2 IS TO BE REPLACED.
!
! MODULES REFERENCED BY SWPTST - NONE
!
!***********************************************************
!
      real x4 , y4 , z4 , dx1 , dx2 , dx3 , dy1 , dy2 , dy3 , dz1 ,     &
     &     dz2 , dz3
!
! LOCAL PARAMETERS -
!
! X4,Y4,Z4 =    COORDINATES OF N4
! DX1,DY1,DZ1 = COORDINATES OF N1 - N4
! DX2,DY2,DZ2 = COORDINATES OF N2 - N4
! DX3,DY3,DZ3 = COORDINATES OF N3 - N4
!
      x4 = x(n4)
      y4 = y(n4)
      z4 = z(n4)
      dx1 = x(n1) - x4
      dx2 = x(n2) - x4
      dx3 = x(n3) - x4
      dy1 = y(n1) - y4
      dy2 = y(n2) - y4
      dy3 = y(n3) - y4
      dz1 = z(n1) - z4
      dz2 = z(n2) - z4
      dz3 = z(n3) - z4
!
! N4 LIES ABOVE THE PLANE OF (N1,N2,N3) IFF N3 LIES ABOVE
!   THE PLANE OF (N2,N1,N4) IFF DET(N3-N4,N2-N4,N1-N4) =
!   (N3-N4,N2-N4 X N1-N4) .GT. 0.
!
      swptst = dx3*(dy2*dz1-dy1*dz2) - dy3*(dx2*dz1-dx1*dz2)            &
     &         + dz3*(dx2*dy1-dx1*dy2)>0.
      end function swptst
!*==TRANS.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine trans(n,rlat,rlon,x,y,z)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - RLAT
!A INPUT  - RLON
!A OUTPUT - X
!A OUTPUT - Y
!A OUTPUT - Z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0001   AA0002
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  COSPHI   I        NN       PHI      THETA
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n
      real rlat(n) , rlon(n) , x(n) , y(n) , z(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE TRANSFORMS SPHERICAL COORDINATES INTO
! CARTESIAN COORDINATES ON THE UNIT SPHERE FOR INPUT TO
! SUBROUTINE TRMESH.  STORAGE FOR X AND Y MAY COINCIDE WITH
! STORAGE FOR RLAT AND RLON IF THE LATTER NEED NOT BE SAVED.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES (POINTS ON THE
!                           UNIT SPHERE) WHOSE COORDINATES
!                           ARE TO BE TRANSFORMED.
!
!                    RLAT - N-VECTOR CONTAINING LATITUDINAL
!                           COORDINATES OF THE NODES IN
!                           RADIANS.
!
!                    RLON - N-VECTOR CONTAINING LONGITUDINAL
!                           COORDINATES OF THE NODES IN
!                           RADIANS.
!
!                   X,Y,Z - VECTORS OF LENGTH .GE. N.
!
! N, RLAT, AND RLON ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - X,Y,Z - CARTESIAN COORDINATES IN THE
!                             RANGE (-1,1).  X(I)**2 +
!                             Y(I)**2 + Z(I)**2 = 1 FOR
!                             I = 1,...,N.
!
! MODULES REFERENCED BY TRANS - NONE
!
! INTRINSIC FUNCTIONS CALLED BY TRANS - COS, SIN
!
!***********************************************************
!
      integer nn , i
      real phi , theta , cosphi
!
! LOCAL PARAMETERS -
!
! NN =     LOCAL COPY OF N
! I =      DO-LOOP INDEX
! PHI =    LATITUDE
! THETA =  LONGITUDE
! COSPHI = COS(PHI)
!
      nn = n
      do i = 1 , nn
         phi = rlat(i)
         theta = rlon(i)
         cosphi = cos(phi)
         x(i) = cosphi*cos(theta)
         y(i) = cosphi*sin(theta)
         z(i) = sin(phi)
      enddo
      end subroutine trans
!*==TRFIND.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine trfind(nst,p,x,y,z,iadj,iend,b1,b2,b3,i1,i2,i3)
!*--********************************************************************
!A INPUT  - NST
!A INPUT  - P
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - IADJ
!A INPUT  - IEND
!A OUTPUT - B1
!A OUTPUT - B2
!A OUTPUT - B3
!A OUTPUT - I1
!A OUTPUT - I2
!A OUTPUT - I3
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   ADNODE   ADNODE   INTRC0   INTRC0   INTRC1   INTRC1
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  B3P1     IND      INDX     LUN      N0       N1       N1S      N2       N2S      N3       N4       NEXT     NF       NL       NRMAX    NRST     PTN1     PTN2     Q        S12      X0       X1       X2       XP       Y0       Y1       Y2
!             YP       Z0       Z1       Z2       ZP
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real det , x0 , x1 , x2 , y0 , y1 , y2 , z0 , z1 , z2
!*** End of declarations inserted by SPAG
      integer nst , iadj(*) , iend(*) , i1 , i2 , i3
      real p(3) , x(*) , y(*) , z(*) , b1 , b2 , b3
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE LOCATES A POINT P RELATIVE TO A TRIANGU-
! LATION OF THE CONVEX HULL OF A SET OF NODES (POINTS ON THE
! UNIT SPHERE).
!
! INPUT PARAMETERS -    NST - INDEX OF NODE AT WHICH TRFIND
!                             BEGINS SEARCH.  SEARCH TIME
!                             DEPENDS ON THE PROXIMITY OF
!                             NST TO P.
!
!                         P - X-, Y-, AND Z-COORDINATES (IN
!                             THAT ORDER) OF THE POINT TO BE
!                             LOCATED.
!
!                     X,Y,Z - VECTORS CONTAINING CARTESIAN
!                             COORDINATES OF THE NODES IN
!                             THE MESH.  (X(I),Y(I),Z(I))
!                             DEFINES NODE I FOR I = 1,...,N
!                             WHERE N .GE. 3.
!
!                      IADJ - SET OF ADJACENCY LISTS OF
!                             NODES IN THE MESH.
!
!                      IEND - POINTERS TO THE ENDS OF
!                             ADJACENCY LISTS IN IADJ FOR
!                             EACH NODE IN THE MESH.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - B1,B2,B3 - UNNORMALIZED BARYCENTRIC
!                                COORDINATES OF THE CENTRAL
!                                PROJECTION OF P ONTO THE
!                                UNDERLYING PLANAR TRIANGLE
!                                IF P IS IN THE CONVEX HULL
!                                OF THE NODES, UNCHANGED IF
!                                I1 = 0.
!
!                     I1,I2,I3 - COUNTERCLOCKWISE-ORDERED
!                                VERTEX INDICES OF A TRI-
!                                ANGLE CONTAINING P IF P IS
!                                AN INTERIOR POINT.  IF P IS
!                                NOT CONTAINED IN THE CONVEX
!                                HULL OF THE NODES, I1 AND
!                                I2 ARE THE RIGHTMOST AND
!                                LEFTMOST NODES WHICH ARE
!                                VISIBLE FROM P, AND I3 = 0.
!                                IF ALL BOUNDARY NODES ARE
!                                VISIBLE FROM P, I1 AND I2
!                                COINCIDE.  IF P AND ALL OF
!                                THE NODES LIE ON A SINGLE
!                                GREAT CIRCLE THEN I1 = I2
!                                = I3 = 0.
!
! MODULES REFERENCED BY TRFIND - NONE
!
! INTRINSIC FUNCTIONS CALLED BY TRFIND - MAX0, ABS
!
!***********************************************************
!
      integer n0 , n1 , n2 , n3 , n4 , indx , ind , nf , nl , n1s ,     &
     &        n2s , next , nrst , nrmax , lun
      real xp , yp , zp , b3p1 , s12 , ptn1 , ptn2 , q(3)
      data nrmax/5/ , lun/6/
!
! LOCAL PARAMETERS -
!
! N0,N1,N2 = NODES IN COUNTERCLOCKWISE ORDER DEFINING A
!              CONE (WITH VERTEX N0) CONTAINING P, OR END-
!              POINTS OF A BOUNDARY EDGE SUCH THAT P RIGHT
!              N1->N2
! N3,N4 =    NODES OPPOSITE N1->N2 AND N2->N1, RESPECTIVELY
! INDX,IND = INDICES FOR IADJ
! NF,NL =    FIRST AND LAST NEIGHBORS OF N0 IN IADJ, OR
!              FIRST (RIGHTMOST) AND LAST (LEFTMOST) NODES
!              VISIBLE FROM P WHEN P IS EXTERIOR TO THE
!              BOUNDARY
! N1S,N2S =  INITIALLY-DETERMINED VALUES OF N1 AND N2 WHEN
!              P IS EXTERIOR
! NEXT =     CANDIDATE FOR I1 OR I2 WHEN P IS EXTERIOR
! NRST =     NUMBER OF RESTARTS WITH NEW N0
! NRMAX =    MAXIMUM ALLOWABLE NUMBER OF RESTARTS BEFORE
!              TERMINATING WITH AN ERROR MESSAGE
! LUN =      LOGICAL UNIT FOR ERROR MESSAGES
! XP,YP,ZP = LOCAL VARIABLES CONTAINING P(1), P(2), AND P(3)
! B3P1 =     B3 + 1 -- B3P1 = 1 IFF B3 = 0 TO WITHIN THE
!              MACHINE PRECISION
! S12 =      SCALAR PRODUCT <N1,N2>
! PTN1 =     SCALAR PRODUCT <P,N1>
! PTN2 =     SCALAR PRODUCT <P,N2>
! Q =        (N2 X N1) X N2  OR  N1 X (N2 X N1) -- USED IN
!              THE BOUNDARY TRAVERSAL WHEN P IS EXTERIOR
! DET =      STATEMENT FUNCTION WHICH COMPUTES A DETERMI-
!              NANT -- DET(X1,...,Z0) .GE. 0 IFF (X0,Y0,Z0)
!              IS IN THE (CLOSED) LEFT HEMISPHERE DEFINED BY
!              THE PLANE CONTAINING (0,0,0), (X1,Y1,Z1),
!              AND (X2,Y2,Z2) WHERE LEFT IS DEFINED RELA-
!              TIVE TO AN OBSERVER AT (X1,Y1,Z1) FACING
!              (X2,Y2,Z2).
      det(x1,y1,z1,x2,y2,z2,x0,y0,z0) = x0*(y1*z2-y2*z1)                &
     &                                  - y0*(x1*z2-x2*z1)              &
     &                                  + z0*(x1*y2-x2*y1)
!
! INITIALIZE VARIABLES
!
      xp = p(1)
      yp = p(2)
      zp = p(3)
      nrst = 0
      n0 = max0(nst,1)
!
! FIND A CONE WITH VERTEX N0 CONTAINING P
!
 100  indx = iend(n0)
      nl = iadj(indx)
      indx = 1
      if ( n0/=1 ) indx = iend(n0-1) + 1
      nf = iadj(indx)
      n1 = nf
      if ( nl/=0 ) then
!
! N0 IS AN INTERIOR NODE.  FIND N1.
!
         do while ( det(x(n0),y(n0),z(n0),x(n1),y(n1),z(n1),xp,yp,zp)   &
     &              <0. )
            indx = indx + 1
            n1 = iadj(indx)
            if ( n1==nl ) then
!
! P IS BETWEEN ARCS N0->N1 AND N0->NF
!
               n2 = nf
!
! P IS CONTAINED IN A CONE DEFINED BY LINE SEGMENTS N0-N1
!   AND N0-N2 WHERE N1 IS ADJACENT TO N2
!
               n3 = n0
               goto 200
            endif
         enddo
      else
!
! N0 IS A BOUNDARY NODE.  SET NL TO THE LAST NONZERO
!   NEIGHBOR OF N0.
!
         ind = iend(n0) - 1
         nl = iadj(ind)
         if ( det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf),xp,yp,zp)<0. )    &
     &        then
!
! P IS TO THE RIGHT OF THE BOUNDARY EDGE N0->NF
!
            n1 = n0
            n2 = nf
            goto 300
         elseif ( det(x(nl),y(nl),z(nl),x(n0),y(n0),z(n0),xp,yp,zp)<0. )&
     &            then
!
! P IS TO THE RIGHT OF THE BOUNDARY EDGE NL->N0
!
            n1 = nl
            n2 = n0
            goto 300
         endif
      endif
      do
!
! P IS TO THE LEFT OF ARCS N0->N1 AND NL->N0.  SET N2 TO THE
!   NEXT NEIGHBOR OF N0 (FOLLOWING N1).
!
         indx = indx + 1
         n2 = iadj(indx)
         if ( det(x(n0),y(n0),z(n0),x(n2),y(n2),z(n2),xp,yp,zp)<0. )    &
     &        then
            n3 = n0
            exit
         else
            n1 = n2
            if ( n1==nl ) then
               if ( det(x(n0),y(n0),z(n0),x(nf),y(nf),z(nf),xp,yp,zp)   &
     &              <0. ) then
                  n2 = nf
                  n3 = n0
                  exit
               else
!
! P IS LEFT OF OR ON ARCS N0->NB FOR ALL NEIGHBORS NB
!   OF N0.  TEST FOR P = +/-N0.
!
                  if ( abs(x(n0)*xp+y(n0)*yp+z(n0)*zp)<1. ) then
!
! ALL POINTS ARE COLLINEAR IFF P IS LEFT OF NB->N0 FOR
!   ALL NEIGHBORS NB OF N0.  SEARCH THE NEIGHBORS OF N0
!   IN REVERSE ORDER.  NOTE -- N1 = NL AND INDX POINTS TO
!   NL.
!
                     do while ( det(x(n1),y(n1),z(n1),x(n0),y(n0),z(n0),&
     &                          xp,yp,zp)>=0. )
                        if ( n1==nf ) then
!
! ALL POINTS ARE COLLINEAR
!
                           i1 = 0
                           i2 = 0
                           i3 = 0
                           return
                        else
                           indx = indx - 1
                           n1 = iadj(indx)
                        endif
                     enddo
                  endif
!
! P IS TO THE RIGHT OF N1->N0, OR P = +/-N0.  SET N0 TO N1
!   AND START OVER.
!
                  n0 = n1
                  nrst = nrst + 1
                  if ( nrst/=nrmax ) goto 100
                  goto 800
               endif
            endif
         endif
      enddo
 200  do
         b3 = det(x(n1),y(n1),z(n1),x(n2),y(n2),z(n2),xp,yp,zp)
         if ( b3>=0. ) then
!
! P IS IN (N1,N2,N3) UNLESS N0, N1, N2, AND P ARE COLLINEAR
!   OR P IS CLOSE TO -N0.
!
            b3p1 = b3 + 1.
            if ( b3p1<=1. ) then
!
! B3 = 0 AND THUS P LIES ON N1->N2. COMPUTE
!   B1 = DET(P,N2 X N1,N2) AND B2 = DET(P,N1,N2 X N1).
!
               b3 = 0.
               s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
               ptn1 = xp*x(n1) + yp*y(n1) + zp*z(n1)
               ptn2 = xp*x(n2) + yp*y(n2) + zp*z(n2)
               b1 = ptn1 - s12*ptn2
               b2 = ptn2 - s12*ptn1
               if ( b1<0. .or. b2<0. ) then
!
! RESTART WITH N0 = 1
!
                  nrst = nrst + 1
                  if ( nrst==nrmax ) goto 800
                  n0 = 1
                  goto 100
               endif
            else
!
! B3 .NE. 0.
!
               b1 = det(x(n2),y(n2),z(n2),x(n3),y(n3),z(n3),xp,yp,zp)
               b2 = det(x(n3),y(n3),z(n3),x(n1),y(n1),z(n1),xp,yp,zp)
               if ( b1<0. .or. b2<0. ) then
!
! RESTART WITH N0 = N3
!
                  nrst = nrst + 1
                  if ( nrst==nrmax ) goto 800
                  n0 = n3
                  goto 100
               endif
            endif
!
! P IS IN (N1,N2,N3)
!
            i1 = n1
            i2 = n2
            i3 = n3
            return
         else
!
! SET N4 TO THE FIRST NEIGHBOR OF N2 FOLLOWING N1
!
            indx = iend(n2)
            if ( iadj(indx)/=n1 ) then
               do
!
! N1 IS NOT THE LAST NEIGHBOR OF N2
!
                  indx = indx - 1
                  if ( iadj(indx)==n1 ) then
                     n4 = iadj(indx+1)
                     if ( n4/=0 ) exit
                     goto 300
                  endif
               enddo
            else
!
! N1 IS THE LAST NEIGHBOR OF N2.
! SET N4 TO THE FIRST NEIGHBOR.
!
               indx = 1
               if ( n2/=1 ) indx = iend(n2-1) + 1
               n4 = iadj(indx)
            endif
!
! DEFINE A NEW ARC N1->N2 WHICH INTERSECTS THE LINE
!   SEGMENT N0-P
!
            if ( det(x(n0),y(n0),z(n0),x(n4),y(n4),z(n4),xp,yp,zp)>=0. )&
     &           then
               n3 = n1
               n1 = n4
            else
               n3 = n2
               n2 = n4
            endif
         endif
      enddo
!
! P RIGHT N1->N2 WHERE N1->N2 IS A BOUNDARY EDGE.
!   SAVE N1 AND N2, AND SET NL = 0 TO INDICATE THAT
!   NL HAS NOT YET BEEN FOUND.
!
 300  n1s = n1
      n2s = n2
      nl = 0
!
!           COUNTERCLOCKWISE BOUNDARY TRAVERSAL
!
 400  indx = 1
      if ( n2/=1 ) indx = iend(n2-1) + 1
      next = iadj(indx)
      if ( det(x(n2),y(n2),z(n2),x(next),y(next),z(next),xp,yp,zp)>=0. )&
     &     then
!
! N2 IS THE RIGHTMOST VISIBLE NODE IF P FORWARD N2->N1
!   OR NEXT FORWARD N2->N1 -- SET Q TO (N2 X N1) X N2
!
         s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
         q(1) = x(n1) - s12*x(n2)
         q(2) = y(n1) - s12*y(n2)
         q(3) = z(n1) - s12*z(n2)
         if ( xp*q(1)+yp*q(2)+zp*q(3)>=0. ) goto 500
         if ( x(next)*q(1)+y(next)*q(2)+z(next)*q(3)>=0. ) goto 500
!
! N1, N2, NEXT, AND P ARE NEARLY COLLINEAR, AND N2 IS
!   THE LEFTMOST VISIBLE NODE
!
         nl = n2
      endif
      n1 = n2
      n2 = next
      if ( n2/=n1s ) goto 400
!
! ALL BOUNDARY NODES ARE VISIBLE
!
      i1 = n1s
      i2 = n1s
      i3 = 0
      return
!
! N2 IS THE RIGHTMOST VISIBLE NODE
!
 500  nf = n2
      if ( nl/=0 ) goto 700
!
! RESTORE INITIAL VALUES OF N1 AND N2, AND BEGIN SEARCH
!   FOR THE LEFTMOST VISIBLE NODE
!
      n2 = n2s
      n1 = n1s
!
!           CLOCKWISE BOUNDARY TRAVERSAL
!
 600  indx = iend(n1) - 1
      next = iadj(indx)
      if ( det(x(next),y(next),z(next),x(n1),y(n1),z(n1),xp,yp,zp)>=0. )&
     &     then
!
! N1 IS THE LEFTMOST VISIBLE NODE IF P OR NEXT IS
!   FORWARD OF N1->N2 -- COMPUTE Q = N1 X (N2 X N1)
!
         s12 = x(n1)*x(n2) + y(n1)*y(n2) + z(n1)*z(n2)
         q(1) = x(n2) - s12*x(n1)
         q(2) = y(n2) - s12*y(n1)
         q(3) = z(n2) - s12*z(n1)
         if ( xp*q(1)+yp*q(2)+zp*q(3)>=0. ) then
!
! N1 IS THE LEFTMOST VISIBLE NODE
!
            nl = n1
            goto 700
         elseif ( x(next)*q(1)+y(next)*q(2)+z(next)*q(3)>=0. ) then
            nl = n1
            goto 700
         else
!
! P, NEXT, N1, AND N2 ARE NEARLY COLLINEAR AND N1 IS THE
!   RIGHTMOST VISIBLE NODE
!
            nf = n1
         endif
      endif
      n2 = n1
      n1 = next
      goto 600
!
! NF AND NL HAVE BEEN FOUND
!
 700  i1 = nf
      i2 = nl
      i3 = 0
      return
!
! MORE THAN NRMAX RESTARTS.  PRINT AN ERROR MESSAGE AND
!   TERMINATE PROCESSING.
!
 800  write (lun,99001) nst , n0 , n1 , n2 , n3 , n4
99001 format ('1','ERROR IN TRFIND -- POSSIBLE INFINITE ',              &
     &        'LOOP DUE TO F.P. ERROR.'//' ','NST,N0,N1,N2,N3,N4 = ',   &
     &        5(i5,', '),i5//' ',                                       &
     &        'TRY REORDERING THE NODES OR CHANGING NST')
      stop
      end subroutine trfind
!*==TRMESH.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine trmesh(n,x,y,z,iadj,iend,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A OUTPUT - IADJ
!A OUTPUT - IEND
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ADNODE
! called by   AA0001   AA0002   AA0003   AA0004
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  IERR     K        NN       X0       X1       X2       Y0       Y1       Y2       Z0       Z1       Z2
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real x0 , x1 , x2 , y0 , y1 , y2 , z0 , z1 , z2
!*** End of declarations inserted by SPAG
      integer n , iadj(*) , iend(n) , ier
      real x(n) , y(n) , z(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CREATES A THIESSEN TRIANGULATION OF THE
! CONVEX HULL OF N ARBITRARILY DISTRIBUTED POINTS (NODES) ON
! THE UNIT SPHERE.  IF THE NODES ARE NOT CONTAINED IN A
! SINGLE HEMISPHERE, THEIR CONVEX HULL IS THE ENTIRE SPHERE,
! AND THERE ARE NO BOUNDARY NODES.  THE TRIANGULATION IS
! OPTIMAL IN THE SENSE THAT IT IS AS NEARLY EQUIANGULAR AS
! POSSIBLE.  TRMESH IS PART OF AN INTERPOLATION PACKAGE
! WHICH ALSO PROVIDES SUBROUTINES TO REORDER THE NODES,
! TRANSFORM COORDINATES, ADD A NEW NODE, DELETE AN ARC, PLOT
! THE MESH, AND PRINT THE DATA STRUCTURE.
!   UNLESS THE NODES ARE ALREADY ORDERED IN SOME REASONABLE
! FASHION, THEY SHOULD BE REORDERED BY SUBROUTINE REORDR FOR
! INCREASED EFFICIENCY BEFORE CALLING TRMESH.  SPHERICAL
! COORDINATES (LATITUDE AND LONGITUDE) MAY BE CONVERTED TO
! CARTESIAN COORDINATES BY SUBROUTINE TRANS.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            N .GE. 3.
!
!                    X,Y,Z - N-VECTORS CONTAINING CARTESIAN
!                            COORDINATES OF DISTINCT NODES.
!                            (X(I),Y(I),Z(I)) DEFINES NODE
!                            I.  X(I)**2 + Y(I)**2 + Z(I)**2
!                            = 1 FOR ALL I.  THE FIRST THREE
!                            NODES MUST NOT BE COLLINEAR
!                            (LIE ON A SINGLE GREAT CIRCLE).
!
!                     IADJ - VECTOR OF LENGTH .GE. 6*(N-1).
!
!                     IEND - VECTOR OF LENGTH .GE. N.
!
! N, X, Y, AND Z ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ - ADJACENCY LISTS OF NEIGHBORS IN
!                            COUNTERCLOCKWISE ORDER.  THE
!                            LIST FOR NODE I+1 FOLLOWS THAT
!                            FOR NODE I. THE VALUE 0 DENOTES
!                            THE BOUNDARY (OR A PSEUDO-NODE
!                            FROM WHICH ALL BOUNDARY NODES
!                            ARE VISIBLE) AND IS ALWAYS THE
!                            LAST NEIGHBOR OF A BOUNDARY
!                            NODE.  IADJ IS UNCHANGED IF IER
!                            .NE. 0.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS (SETS OF
!                            NEIGHBORS) IN IADJ.  THE
!                            NEIGHBORS OF NODE 1 BEGIN IN
!                            IADJ(1).  FOR K .GT. 1, THE
!                            NEIGHBORS OF NODE K BEGIN IN
!                            IADJ(IEND(K-1)+1) AND K HAS
!                            IEND(K) - IEND(K-1) NEIGHBORS
!                            INCLUDING (POSSIBLY) THE
!                            BOUNDARY.  IADJ(IEND(K)) .EQ. 0
!                            IFF NODE K IS ON THE BOUNDARY.
!                            IEND IS UNCHANGED IF IER .NE.
!                            0.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE
!                                    ENCOUNTERED.
!                            IER = 1 IF N .LT. 3.
!                            IER = 2 IF THE FIRST THREE
!                                    NODES ARE COLLINEAR.
!
! MODULES REFERENCED BY TRMESH - ADNODE, TRFIND, INTADD,
!                                BDYADD, COVSPH, SHIFTD,
!                                INDEX, SWPTST, SWAP
!
!***********************************************************
!
      integer nn , k , ierr
      logical left
!
! LOCAL PARAMETERS -
!
! NN =   LOCAL COPY OF N
! K =    INDEX OF NODE TO BE ADDED, DO-LOOP INDEX
! IERR = ERROR FLAG FOR CALL TO ADNODE - NOT CHECKED
! LEFT = STATEMENT FUNCTION -- TRUE IFF (X0,Y0,Z0) IS IN
!          THE (CLOSED) LEFT HEMISPHERE DEFINED BY THE
!          PLANE CONTAINING (0,0,0), (X1,Y1,Z1), AND
!          (X2,Y2,Z2) WHERE LEFT IS DEFINED RELATIVE TO
!          AN OBSERVER AT (X1,Y1,Z1) FACING (X2,Y2,Z2)
!
      left(x1,y1,z1,x2,y2,z2,x0,y0,z0) = x0*(y1*z2-y2*z1)               &
     &   - y0*(x1*z2-x2*z1) + z0*(x1*y2-x2*y1)>=0.
      nn = n
      ier = 1
      if ( nn<3 ) return
      ier = 0
      if ( .not.left(x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3)) )    &
     &     then
!
! FIRST TRIANGLE IS (3,2,1) = (2,1,3) = (1,3,2).  INITIALIZE
!   IADJ.
!
         iadj(1) = 3
         iadj(2) = 2
         iadj(3) = 0
         iadj(4) = 1
         iadj(5) = 3
         iadj(6) = 0
         iadj(7) = 2
         iadj(8) = 1
         iadj(9) = 0
      elseif ( .not.left(x(2),y(2),z(2),x(1),y(1),z(1),x(3),y(3),z(3)) )&
     &         then
!
! FIRST TRIANGLE IS (1,2,3), I.E. 3 STRICTLY LEFT 1-2,
!   I.E. 3 LIES IN THE LEFT HEMISPHERE DEFINED BY ARC 1-2
!
         iadj(1) = 2
         iadj(2) = 3
         iadj(3) = 0
         iadj(4) = 3
         iadj(5) = 1
         iadj(6) = 0
         iadj(7) = 1
         iadj(8) = 2
         iadj(9) = 0
      else
!
! FIRST 3 NODES ARE COLLINEAR
!
         ier = 2
         return
      endif
!
! INITIALIZE IEND
!
      iend(1) = 3
      iend(2) = 6
      iend(3) = 9
      if ( nn==3 ) return
!
! ADD NODES 4,...,N TO THE TRIANGULATION
!
      do k = 4 , nn
         call adnode(k,x,y,z,iadj,iend,ierr)
      enddo
      end subroutine trmesh
!*==TRPLOT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine trplot(n,x,y,z,iadj,iend,elat,elon,ititle,nc,numbr,ier)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - IADJ
!A OUTPUT - IEND
!A INPUT  - ELAT
!A INPUT  - ELON
!A UNUSED - ITITLE
!A INPUT  - NC
!A UNUSED - NUMBR
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       CIRCLE   PRJCT
! called by   ** NOTHING **
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  BDRY     CX       CY       ER       EX       EY       EZ       I        IERR     INDF     INDL     ISW      N0       N1       NEW      NN       NPC      ST0      TOL      VX       VY       VZ       X0       X1       Y0       Y1
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iadj(*) , iend(n) , ititle(*) , nc , numbr , ier
      real x(n) , y(n) , z(n) , elat , elon
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE PLOTS THE ARCS OF A TRIANGULATION OF A
! SET OF NODES ON THE UNIT SPHERE.  VISIBLE NODES ARE PRO-
! JECTED ONTO A PLANE WHICH CONTAINS THE ORIGIN AND HAS
! NORMAL GIVEN BY A USER-SPECIFIED EYE-POSITION.  PROJEC-
! TIONS OF ADJACENT (VISIBLE) NODES ARE CONNECTED BY LINE
! SEGMENTS AND A UNIT CIRCLE (THE PROJECTION OF THE SPHERE
! BOUNDARY) IS DRAWN AROUND THE (PROJECTED) TRIANGULATION.
! A SPECIAL SYMBOL IS PLOTTED IN THE POSITION OF A NODE
! WHOSE NEIGHBORS ARE NOT ALL VISIBLE FROM THE SPECIFIED
! EYE-POSITION.  CARDS WITH C* IN THE FIRST TWO COLUMNS
! MUST BE REPLACED WITH CALLS TO USER-SUPPLIED GRAPHICS
! SUBROUTINES IN ORDER TO MAKE USE OF THIS ROUTINE.
!
! INPUT PARAMETERS -    N - NUMBER OF NODES.  N .GE. 3.
!
!                   X,Y,Z - N-VECTORS CONTAINING THE
!                           CARTESIAN COORDINATES OF THE
!                           NODES.
!
!               IADJ,IEND - TRIANGULATION DATA STRUCTURE
!                           (SEE SUBROUTINE TRMESH).
!
!               ELAT,ELON - LATITUDE AND LONGITUDE OF THE
!                           ENDPOINT OF THE UNIT NORMAL TO
!                           THE PROJECTION PLANE.  THE
!                           VECTOR FROM THE ORIGIN (0,0,0)
!                           TO THE EYE-POSITION IS TAKEN TO
!                           BE A LARGE MULTIPLE OF THIS UNIT
!                           NORMAL VECTOR.  THUS ELAT AND
!                           ELON SPECIFY THE POINT AT THE
!                           CENTER OF THE PLOT.  -PI/2 .LE.
!                           ELAT .LE. PI/2 AND -PI .LE. ELON
!                           .LE. PI.
!
!                  ITITLE - INTEGER ARRAY CONTAINING A LINE
!                           OF TEXT TO BE CENTERED ABOVE THE
!                           PLOT IF NC .GT. 0.  ITITLE MUST
!                           BE INITIALIZED WITH HOLLERITH
!                           CONSTANTS OR READ WITH AN A-FOR-
!                           MAT.  ITS DIMENSION DEPENDS ON
!                           NC AND THE NUMBER OF CHARACTERS
!                           STORED IN A COMPUTER WORD.
!
!                      NC - NUMBER OF CHARACTERS IN ITITLE.
!                           0 .LE. NC .LE. 40.  NO TITLE IS
!                           DRAWN IF NC = 0.
!
!                   NUMBR - OPTION INDICATOR.  IF NUMBR .NE.
!                           0, THE NODAL INDICES ARE PLOTTED
!                           NEXT TO THE NODES.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER  - IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE
!                                   ENCOUNTERED.
!                           IER = 1 IF N OR NC IS OUT OF
!                                   RANGE.
!
!   THE VALUES IN THE DATA STATEMENT BELOW MAY BE ALTERED
! IN ORDER TO MODIFY VARIOUS PLOTTING OPTIONS.
!
! MODULES REFERENCED BY TRPLOT - CIRCLE, PRJCT
!
! INTRINSIC FUNCTIONS CALLED BY TRPLOT - COS, SIN, ATAN,
!                                        ABS, IABS
!
!***********************************************************
!
      integer npc , nn , isw , n0 , indf , indl , i , ierr , n1 , new
      real er , tol , cx(481) , cy(481) , ex , ey , ez , vx , vy , vz , &
     &     x0 , y0 , x1 , y1
      logical st0 , bdry
      data npc , er , tol/481 , 50. , .004/
!
! LOCAL PARAMETERS -
!
! NPC =      DIMENSION OF CX AND CY, AND NUMBER OF POINTS
!              (PLUS 1) USED TO PLOT THE UNIT CIRCLE --
!              DETERMINES HOW CLOSELY THE CIRCLE IS APPROX-
!              IMATED.  NOTE THAT IF NPC IS INCREASED, THE
!              DIMENSION STATEMENTS MUST ALSO BE INCREASED.
! ER =       RADIAL COORDINATE OF THE EYE-POSITION (DISTANCE
!              FROM THE ORIGIN)
! TOL =      MINIMUM ALLOWABLE ARC-LENGTH BETWEEN THE EYE-
!              POSITION AND (VX,VY,VZ)
! NN =       LOCAL COPY OF N
! ISW =      SWITCH USED BY SUBROUTINE PRJCT
! N0 =       NODE WHOSE INCIDENT ARCS ARE TO BE PLOTTED
! INDF =     IADJ INDEX OF THE FIRST NEIGHBOR OF N0
! INDL =     IADJ INDEX OF THE LAST NEIGHBOR OF N0
! I =        IADJ INDEX IN THE RANGE INDF,...,INDL
! IERR =     ERROR FLAG FOR CALLS TO PRJCT -- NOT CHECKED
! N1 =       NEIGHBOR OF N0
! NEW =      NEIGHBOR OF N0 AND CANDIDATE FOR NEXT VALUE OF
!              N0
! CX,CY =    ARRAYS CONTAINING THE COORDINATES OF POINTS ON
!              THE UNIT CIRCLE
! EX,EY,EZ = CARTESIAN COORDINATES OF THE EYE-POSITION
! VX,VY,VZ = CARTESIAN COORDINATES OF A POINT WHOSE PROJECT-
!              ION DETERMINES THE POSITIVE Y-AXIS FOR THE
!              PLOT
! X0,Y0 =    PROJECTED COORDINATES OF N0
! X1,Y1 =    PROJECTED COORDINATES OF N1
! ST0 =      SWITCH USED TO ALTERNATE DIRECTION OF PEN
!              MOVEMENT
! BDRY =     SWITCH DETERMINING WHETHER ALL NEIGHBORS OF N0
!              ARE VISIBLE
!
      nn = n
!
! CHECK FOR INVALID PARAMETERS
!
      ier = 1
      if ( nn<3 .or. nc<0 .or. nc>40 ) return
      ier = 0
!
! COMMANDS WHICH PERFORM THE FOLLOWING FUNCTIONS SHOULD BE
! INSERTED HERE --
!
!*  INITIALIZE THE PLOTTING ENVIRONMENT IF NECESSARY,
!*  SET DIMENSIONS OF THE PLOTTER SPACE,
!*  ESTABLISH A LINEAR MAPPING FROM THE DATA SPACE (-1,1) X
!*    (-1,1) TO THE PLOTTER SPACE, AND
!*  PLOT THE TITLE IF NC .NE. 0.
!
! PLOT THE UNIT CIRCLE
!
      call circle(npc,cx,cy)
!*    CALL CURVE(CX,CY,NPC)
!
! SET PARAMETERS FOR PRJCT.  THE Y-AXIS WILL BE THE PRO-
!   JECTION OF THE NORTH POLE UNLESS THAT IS TOO CLOSE TO
!   THE EYE-POSITION.
!
      ex = er*cos(elon)*cos(elat)
      ey = er*sin(elon)*cos(elat)
      ez = er*sin(elat)
      vx = 0.
      vy = 0.
      vz = 1.
      if ( 2.*atan(1.)-abs(elat)<tol ) then
         vx = 1.
         if ( ez>0. ) vx = -1.
         vz = 0.
      endif
      isw = 1
!
! INITIALIZE FOR LOOP ON NODES.  EACH VISIBLE NODE N0 IS
!   CONNECTED TO ITS (VISIBLE) NEIGHBORS WHICH HAVE LARGER
!   INDICES.  N0 IS THEN MARKED BY MAKING THE CORRESPONDING
!   IEND ENTRY NEGATIVE, AND THE SEARCH FOR THE NEXT
!   UNMARKED NODE BEGINS WITH THE NEIGHBORS OF N0.
!
      n0 = 1
      indf = 1
!
! TOP OF LOOP -- SET INDL AND TEST FOR N0 VISIBLE
!
 100  indl = iend(n0)
      if ( ex*x(n0)+ey*y(n0)+ez*z(n0)>=0. ) then
!
! INITIALIZE ST0 AND BDRY, COMPUTE X0,Y0, AND INITIALIZE I
!   FOR NEIGHBOR LOOP
!
         st0 = .true.
         bdry = .false.
         call prjct(x(n0),y(n0),z(n0),0.,0.,0.,ex,ey,ez,vx,vy,vz,isw,x0,&
     &              y0,ierr)
         i = indl
         if ( iadj(i)==0 ) i = i - 1
      else
         iend(n0) = -indl
         i = indf
         goto 300
      endif
!
! LOOP ON NEIGHBORS OF N0 IN REVERSE ORDER
!
 200  n1 = iadj(i)
      if ( ex*x(n1)+ey*y(n1)+ez*z(n1)<0. ) then
!
! N1 NOT VISIBLE -- SET BDRY AND MARK N1
!
         bdry = .true.
         iend(n1) = -iabs(iend(n1))
!
! N1 VISIBLE
!
      elseif ( n1>=n0 ) then
!
! N1 IS VISIBLE AND HAS LARGER INDEX THAN N0 -- COMPUTE ITS
!   PROJECTION
!
         call prjct(x(n1),y(n1),z(n1),0.,0.,0.,ex,ey,ez,vx,vy,vz,isw,x1,&
     &              y1,ierr)
!
! CONNECT N0 AND N1 -- THE DIRECTION OF PEN MOVEMENT ALTER-
!   NATES BETWEEN AWAY FROM N0 AND TOWARD N0 FOR REDUCED
!   PEN-UP TIME
!
!*    IF (ST0) CALL LINE(X0,Y0,X1,Y1)
!*    IF (.NOT. ST0) CALL LINE(X1,Y1,X0,Y0)
         st0 = .not.st0
      endif
!
! TEST FOR TERMINATION OF NEIGHBOR LOOP
!
      if ( i==indf ) then
!
! MARK N0 AS HAVING BEEN PROCESSED, PLOT A SYMBOL IF ANY OF
!   ITS NEIGHBORS ARE NOT VISIBLE, AND PLOT THE INDEX IF
!   NUMBR .NE. 0.
!
         iend(n0) = -indl
      else
         i = i - 1
         goto 200
      endif
 300  do
!*    IF (BDRY) CALL POINT(X0,Y0)
!*    IF (NUMBR .NE. 0) CALL RLINT(N0,X0,Y0)
!
! SEARCH THE NEIGHBORS OF N0 FOR AN UNMARKED NODE STARTING
!   WITH IADJ(I) = IADJ(INDF)
!
         new = iadj(i)
         if ( new==0 ) exit
         if ( iend(new)<0 ) then
!
! TEST FOR TERMINATION
!
            if ( i==indl ) exit
            i = i + 1
         else
            n0 = new
            indf = iabs(iend(n0-1)) + 1
            goto 100
         endif
      enddo
!
! ALL NEIGHBORS OF N0 ARE MARKED.  SEARCH IEND FOR AN
!   UNMARKED NODE
!
      do n0 = 2 , nn
         if ( iend(n0)>=0 ) then
            indf = -iend(n0-1) + 1
            goto 100
         endif
      enddo
!
! ALL NODES HAVE BEEN MARKED -- RESTORE IEND
!
      do n0 = 1 , nn
         iend(n0) = -iend(n0)
      enddo
!
! TERMINATE PLOTTING -- MOVE TO A NEW FRAME
!
!*    CALL FRAME
      end subroutine trplot
!*==TRPRNT.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine trprnt(n,lunit,x,y,z,iadj,iend,iflag)
!*--********************************************************************
!A INPUT  - N
!A INPUT  - LUNIT
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - Z
!A INPUT  - IADJ
!A INPUT  - IEND
!A INPUT  - IFLAG
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0002   AA0004
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        INC      INDF     INDL     LUN      NA       NB       NL       NLMAX    NMAX     NN       NODE     NT
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , lunit , iadj(*) , iend(n) , iflag
      real x(n) , y(n) , z(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF POINTS ON THE UNIT
! SPHERE, THIS ROUTINE PRINTS THE ADJACENCY LISTS AND,
! OPTIONALLY, THE NODAL COORDINATES (EITHER LATITUDE AND
! LONGITUDE OR CARTESIAN COORDINATES).  THE NUMBERS OF
! BOUNDARY NODES, TRIANGLES, AND ARCS ARE ALSO PRINTED.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            3 .LE. N .LE. 9999.
!
!                    LUNIT - LOGICAL UNIT FOR OUTPUT.  1
!                            .LE. LUNIT .LE. 99.  OUTPUT IS
!                            PRINTED ON UNIT 6 IF LUNIT IS
!                            OUT OF RANGE.
!
!                    X,Y,Z - VECTORS OF COORDINATES OF THE
!                            NODES IN THE MESH.  THESE MAY
!                            BE DUMMY PARAMETERS IF IFLAG IS
!                            NOT 0 OR 1, AND Z IS NOT USED
!                            IF IFLAG .NE. 0.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF X, Y, AND Z ARE TO
!                                      BE PRINTED (TO 6 DEC-
!                                      IMAL PLACES).
!                            IFLAG = 1 IF ONLY X AND Y (AS-
!                                      SUMED TO CONTAIN LON-
!                                      GITUDE AND LATITUDE,
!                                      RESPECTIVELY) ARE TO
!                                      BE PRINTED.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! NONE OF THE PARAMETERS ARE ALTERED BY THIS ROUTINE.
!
! MODULES REFERENCED BY TRPRNT - NONE
!
!***********************************************************
!
      integer nn , nmax , lun , node , indf , indl , nl , nlmax , inc , &
     &        i , nb , nt , na
      data nmax/9999/ , nlmax/60/
!
! LOCAL PARAMETERS -
!
! NN =        LOCAL COPY OF N
! NMAX =      UPPER BOUND ON N
! LUN =       LOCAL COPY OF LUNIT
! NODE =      INDEX OF A NODE
! INDF,INDL = IADJ INDICES OF THE FIRST AND LAST NEIGHBORS
!               OF NODE
! NL =        NUMBER OF LINES PRINTED ON A PAGE
! NLMAX =     MAXIMUM NUMBER OF PRINT LINES PER PAGE EXCEPT
!               FOR THE LAST PAGE WHICH MAY HAVE 2 ADDITION-
!               AL LINES
! INC =       INCREMENT FOR NL
! I =         IADJ INDEX FOR IMPLIED DO-LOOP
! NB =        NUMBER OF BOUNDARY NODES
! NT =        NUMBER OF TRIANGLES
! NA =        NUMBER OF ARCS (UNDIRECTED EDGES)
!
      nn = n
      lun = lunit
      if ( lun<1 .or. lun>99 ) lun = 6
!
! PRINT HEADING AND INITIALIZE NL
!
      write (lun,99001) nn
!
! PRINT FORMATS
!
99001 format ('1',26x,'ADJACENCY SETS,    N = ',i5//)
      if ( nn<3 .or. nn>nmax ) then
!
! N IS OUT OF RANGE
!
         write (lun,99002)
99002    format (' ',10x,'*** N IS OUT OF RANGE ***')
         return
      else
         nl = 6
         if ( iflag==0 ) then
!
! PRINT X, Y, Z, AND IADJ
!
            write (lun,99003)
99003       format (' ','NODE',5x,'X(NODE)',8x,'Y(NODE)',8x,'Z(NODE)',  &
     &              12x,'NEIGHBORS OF NODE'//)
            nb = 0
            indf = 1
            do node = 1 , nn
               indl = iend(node)
               if ( iadj(indl)==0 ) nb = nb + 1
               inc = (indl-indf)/5 + 2
               nl = nl + inc
               if ( nl>nlmax ) write (lun,99011)
               if ( nl>nlmax ) nl = inc
               write (lun,99004) node , x(node) , y(node) , z(node) ,   &
     &                           (iadj(i),i=indf,indl)
99004          format (' ',i4,3E15.6,5x,5I5/(' ',54x,5I5))
               if ( indl-indf/=4 ) write (lun,99010)
               indf = indl + 1
            enddo
         elseif ( iflag==1 ) then
!
! PRINT X (LONGITUDE), Y (LATITUDE), AND IADJ
!
            write (lun,99005)
99005       format (' ','NODE',4x,'LONGITUDE',7x,'LATITUDE',19x,        &
     &              'NEIGHBORS OF NODE'//)
            nb = 0
            indf = 1
            do node = 1 , nn
               indl = iend(node)
               if ( iadj(indl)==0 ) nb = nb + 1
               inc = (indl-indf)/8 + 2
               nl = nl + inc
               if ( nl>nlmax ) write (lun,99011)
               if ( nl>nlmax ) nl = inc
               write (lun,99006) node , x(node) , y(node) ,             &
     &                           (iadj(i),i=indf,indl)
99006          format (' ',i4,2E15.6,5x,8I5/(' ',39x,8I5))
               if ( indl-indf/=7 ) write (lun,99010)
               indf = indl + 1
            enddo
         else
!
! PRINT IADJ ONLY
!
            write (lun,99007)
99007       format (' ','NODE',32x,'NEIGHBORS OF NODE'//)
            nb = 0
            indf = 1
            do node = 1 , nn
               indl = iend(node)
               if ( iadj(indl)==0 ) nb = nb + 1
               inc = (indl-indf)/14 + 2
               nl = nl + inc
               if ( nl>nlmax ) write (lun,99011)
               if ( nl>nlmax ) nl = inc
               write (lun,99008) node , (iadj(i),i=indf,indl)
99008          format (' ',i4,5x,14I5/(' ',9x,14I5))
               if ( indl-indf/=13 ) write (lun,99010)
               indf = indl + 1
            enddo
         endif
      endif
!
! PRINT NB, NA, AND NT
!
      nt = 2*nn - nb - 2
      if ( nb==0 ) nt = nt - 2
      na = nt + nn - 1
      if ( nb==0 ) na = na - 1
      write (lun,99009) nb , na , nt
99009 format (/' ','NB = ',i4,' BOUNDARY NODES',10x,'NA = ',i5,' ARCS', &
     &        10x,'NT = ',i5,' TRIANGLES')
      return
99010 format (' ')
99011 format ('1')
      end subroutine trprnt
!*==UNIF.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine unif(n,x,y,z,w,iadj,iend,nrow,ni,nj,plat,plon,iflag,   &
     &                grad,ww,ier)
!*--********************************************************************
!A INPUT  - N
!A PASSED - X
!A PASSED - Y
!A PASSED - Z
!A PASSED - W
!A PASSED - IADJ
!A PASSED - IEND
!A INPUT  - NROW
!A INPUT  - NI
!A INPUT  - NJ
!A PASSED - PLAT
!A PASSED - PLON
!A INPUT  - IFLAG
!A PASSED - GRAD
!A PASSED - WW
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       GRADL    INTRC1
! called by   AA0001   AA0002   AA0003   AA0004
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  I        IERR     IFL      IST      J        NEX      NN       NST      NX       NY
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , iadj(*) , iend(n) , nrow , ni , nj , iflag , ier
      real x(n) , y(n) , z(n) , w(n) , plat(ni) , plon(nj) , grad(3,n) ,&
     &     ww(nrow,nj)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF A SET OF NODES ON THE
! UNIT SPHERE AND CORRESPONDING DATA VALUES, THIS ROUTINE
! INTERPOLATES THE DATA VALUES TO A UNIFORM GRID FOR SUCH
! APPLICATIONS AS CONTOURING.  THE INTERPOLANT IS ONCE CON-
! TINUOUSLY DIFFERENTIABLE.  EXTRAPOLATION IS PERFORMED AT
! GRID POINTS EXTERIOR TO THE TRIANGULATION WHEN THE NODES
! DO NOT COVER THE ENTIRE SPHERE.
!
! INPUT PARAMETERS - N - NUMBER OF NODES.  N .GE. 3 AND
!                        N .GE. 7 IF IFLAG .NE. 1.
!
!                X,Y,Z - VECTORS CONTAINING CARTESIAN COORD-
!                        INATES OF THE NODES.
!
!                    W - VECTOR CONTAINING DATA VALUES.
!                        W(I) IS ASSOCIATED WITH (X(I),
!                        Y(I),Z(I)).
!
!            IADJ,IEND - TRIANGULATION DATA STRUCTURE - MAY
!                        BE CREATED BY TRMESH.
!
!                 NROW - NUMBER OF ROWS IN THE DIMENSION
!                        STATEMENT OF WW.
!
!                NI,NJ - NUMBER OF ROWS AND COLUMNS IN THE
!                        UNIFORM GRID.  1 .LE. NI .LE. NROW,
!                        1 .LE. NJ.
!
!            PLAT,PLON - VECTORS OF LENGTH NI AND NJ, RE-
!                        SPECTIVELY, CONTAINING THE LATI-
!                        TUDES AND LONGITUDES OF THE GRID
!                        LINES.
!
!                IFLAG - OPTION INDICATOR
!                        IFLAG = 0 IF GRADIENT ESTIMATES AT
!                                  THE VERTICES OF A TRIAN-
!                                  GLE ARE TO BE RECOMPUTED
!                                  FOR EACH GRID POINT IN
!                                  THE TRIANGLE AND NOT
!                                  SAVED.
!                        IFLAG = 1 IF GRADIENT ESTIMATES ARE
!                                  INPUT IN GRAD.
!                        IFLAG = 2 IF GRADIENT ESTIMATES ARE
!                                  TO BE COMPUTED ONCE FOR
!                                  EACH NODE (BY GRADL) AND
!                                  SAVED IN GRAD.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!                 GRAD - NOT USED IF IFLAG = 0, 3 BY N ARRAY
!                        WHOSE COLUMNS CONTAIN THE X, Y, AND
!                        Z COMPONENTS (IN THAT ORDER) OF THE
!                        GRADIENTS AT THE NODES IF IFLAG =
!                        1, OR ARRAY OF SUFFICIENT SIZE IF
!                        IFLAG = 2.
!
! GRADIENT ESTIMATES MAY BE COMPUTED BY GRADL OR GRADG IF
!   IFLAG = 1.
!
!                   WW - NROW BY NCOL ARRAY WITH NROW .GE.
!                        NI AND NCOL .GE. NJ.
!
! OUTPUT PARAMETERS - GRAD - ARRAY CONTAINING ESTIMATED
!                            GRADIENTS AS DEFINED ABOVE IF
!                            IFLAG = 2 AND IER .GE. 0.  GRAD
!                            IS NOT ALTERED IF IFLAG .NE. 2.
!
!                       WW - INTERPOLATED VALUES AT THE GRID
!                            POINTS IF IER .GE. 0.  WW(I,J)
!                            = F(PLAT(I),PLON(J)) FOR I =
!                            1,...,NI AND J = 1,...,NJ.
!
!                      IER - ERROR INDICATOR
!                            IER = K IF IF NO ERRORS WERE
!                                    ENCOUNTERED AND K GRID
!                                    POINTS REQUIRED EX-
!                                    TRAPOLATION FOR K .GE.
!                                    0.
!                            IER = -1 IF N, NI, NJ, OR IFLAG
!                                     IS OUT OF RANGE.
!                            IER = -2 IF THE NODES ARE COL-
!                                     LINEAR.
!                            IER = -3 IF EXTRAPOLATION FAIL-
!                                     ED DUE TO THE UNIFORM
!                                     GRID EXTENDING TOO FAR
!                                     BEYOND THE TRIANGULA-
!                                     TION BOUNDARY.
!
! MODULES REFERENCED BY UNIF - INTRC1, TRFIND, WVAL, ARCINT,
!                              ARCLEN,
!           (AND OPTIONALLY)   GRADL, GETNP, CONSTR, APLYR,
!                              SETUP, GIVENS, ROTATE, APLYRT
!
!***********************************************************
!
      integer nst , ist , nn , nx , ny , ifl , i , j , ierr , nex
      data nst/1/
!
! LOCAL PARAMETERS -
!
! NST =   INITIAL VALUE FOR IST
! IST =   PARAMETER FOR INTRC1
! NN =    LOCAL COPY OF N
! NX,NY = LOCAL COPIES OF NI AND NJ
! IFL =   LOCAL COPY OF IFLAG
! I,J =   DO-LOOP INDICES
! IERR =  ERROR FLAG FOR CALLS TO GRADL AND INTRC1
! NEX =   NUMBER OF GRID POINTS EXTERIOR TO THE TRIANGULA-
!           TION BOUNDARY (NUMBER OF EXTRAPOLATED VALUES)
!
      nn = n
      nx = ni
      ny = nj
      ifl = iflag
      if ( nx<1 .or. nx>nrow .or. ny<1 .or. ifl<0 .or. ifl>2 ) then
!
! NI, NJ, OR IFLAG IS OUT OF RANGE
!
         ier = -1
         return
      else
         ist = nst
         if ( ifl==2 ) then
!
! COMPUTE GRADIENT ESTIMATES AT THE NODES
!
            do i = 1 , nn
               call gradl(nn,i,x,y,z,w,iadj,iend,grad(1,i),ierr)
               if ( ierr<0 ) goto 100
            enddo
            ifl = 1
         endif
!
! COMPUTE UNIFORM GRID POINTS AND INTERPOLATED VALUES
!
         nex = 0
         do j = 1 , ny
            do i = 1 , nx
               call intrc1(nn,plat(i),plon(j),x,y,z,w,iadj,iend,ifl,    &
     &                     grad,ist,ww(i,j),ierr)
               if ( ierr<0 ) goto 100
               nex = nex + ierr
            enddo
         enddo
         ier = nex
         return
      endif
!
! ERROR IN GRADL OR INTRC1
!
 100  ier = ierr
      end subroutine unif
!*==WVAL.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine wval(b1,b2,b3,v1,v2,v3,w1,w2,w3,g1,g2,g3,iflag,pw,pg)
!*--********************************************************************
!A INPUT  - B1
!A INPUT  - B2
!A INPUT  - B3
!A INPUT  - V1
!A INPUT  - V2
!A INPUT  - V3
!A INPUT  - W1
!A INPUT  - W2
!A INPUT  - W3
!A INPUT  - G1
!A INPUT  - G2
!A INPUT  - G3
!A INPUT  - IFLAG
!A OUTPUT - PW
!A OUTPUT - PG
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ARCINT
! called by   INTRC1
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  C1       C2       C3       DUM      G        I        Q1       Q2       Q3       SUM      U1       U1N      U2       U2N      U3       U3N      VAL      W
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer iflag
      real b1 , b2 , b3 , v1(3) , v2(3) , v3(3) , w1 , w2 , w3 , g1(3) ,&
     &     g2(3) , g3(3) , pw , pg(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN DATA VALUES AND GRADIENTS AT THE THREE VERTICES OF
! A SPHERICAL TRIANGLE CONTAINING A POINT P, THIS ROUTINE
! COMPUTES THE VALUE AND, OPTIONALLY, THE GRADIENT OF F AT P
! WHERE F INTERPOLATES THE VERTEX DATA.  ALONG THE TRIANGLE
! EDGES, THE INTERPOLATORY FUNCTION F IS THE HERMITE CUBIC
! (WITH RESPECT TO ARC-LENGTH) INTERPOLANT OF THE VALUES AND
! TANGENTIAL GRADIENT COMPONENTS AT THE ENDPOINTS, AND THE
! DERIVATIVE NORMAL TO THE ARC VARIES LINEARLY WITH RESPECT
! TO ARC-LENGTH BETWEEN THE NORMAL GRADIENT COMPONENTS AT
! THE ENDPOINTS.  THUS THE METHOD YIELDS C-1 CONTINUITY WHEN
! USED TO INTERPOLATE OVER A TRIANGULATION.  THE INTERPOLANT
! USES A FIRST-ORDER BLENDING METHOD ON THE UNDERLYING
! PLANAR TRIANGLE.
!
! INPUT PARAMETERS - B1,B2,B3 - BARYCENTRIC COORDINATES OF
!                               PP WITH RESPECT TO THE
!                               (PLANAR) UNDERLYING TRIANGLE
!                               (V1,V2,V3) WHERE PP IS THE
!                               CENTRAL PROJECTION OF P ONTO
!                               THIS TRIANGLE.
!
!                    V1,V2,V3 - CARTESIAN COORDINATES OF THE
!                               VERTICES OF A SPHERICAL TRI-
!                               ANGLE CONTAINING P.  V3 LEFT
!                               V1->V2.
!
!                    W1,W2,W3 - DATA VALUES ASSOCIATED WITH
!                               THE VERTICES.
!
!                    G1,G2,G3 - GRADIENTS ASSOCIATED WITH
!                               THE VERTICES.  GI IS ORTHOG-
!                               ONAL TO VI FOR I = 1,2,3.
!
!                       IFLAG - OPTION INDICATOR
!                               IFLAG = 0 IF ONLY PW IS TO
!                                         BE COMPUTED.
!                               IFLAG = 1 IF BOTH PW AND PG
!                                         ARE DESIRED.
!
!                          PG - VECTOR OF LENGTH 3 IF IFLAG
!                               = 1, NOT USED OTHERWISE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -      PW - INTERPOLATED VALUE AT P.
!
!                          PG - INTERPOLATED GRADIENT AT P
!                               (ORTHOGONAL TO P) IF IFLAG
!                               = 1.
!
! EACH VECTOR V ABOVE CONTAINS X-, Y-, AND Z-COMPONENTS IN
!   V(1), V(2), AND V(3), RESPECTIVELY.
!
! MODULES REFERENCED BY WVAL - ARCINT, ARCLEN
!
! INTRINSIC FUNCTION CALLED BY WVAL - SQRT
!
!***********************************************************
!
      integer i
      real c1 , c2 , c3 , sum , u1(3) , u2(3) , u3(3) , u1n , u2n ,     &
     &     u3n , q1(3) , q2(3) , q3(3) , val , w , g(3) , dum
!
! LOCAL PARAMETERS -
!
! I =           DO-LOOP INDEX
! C1,C2,C3 =    COEFFICIENTS (WEIGHT FUNCTIONS) OF PARTIAL
!                 INTERPOLANTS.  C1 = 1 ON THE EDGE OPPOSITE
!                 V1 AND C1 = 0 ON THE OTHER EDGES.  SIMI-
!                 LARLY FOR C2 AND C3.  C1+C2+C3 = 1.
! SUM =         QUANTITY USED TO NORMALIZE C1, C2, AND C3
! U1,U2,U3 =    POINTS ON THE BOUNDARY OF THE PLANAR TRIAN-
!                 GLE AND LYING ON THE LINES CONTAINING PP
!                 AND THE VERTICES.  U1 IS OPPOSITE V1, ETC.
! U1N,U2N,U3N = QUANTITIES USED TO COMPUTE Q1, Q2, AND Q3
!                 (MAGNITUDES OF U1, U2, AND U3)
! Q1,Q2,Q3 =    CENTRAL PROJECTIONS OF U1, U2, AND U3 ONTO
!                 THE SPHERE AND THUS LYING ON AN ARC OF THE
!                 SPHERICAL TRIANGLE
! VAL =         LOCAL VARIABLE USED TO ACCUMULATE THE CON-
!                 TRIBUTIONS TO PW
! W,G =         VALUE AND GRADIENT AT Q1, Q2, OR Q3 OBTAINED
!                 BY INTERPOLATION ALONG ONE OF THE ARCS OF
!                 THE SPHERICAL TRIANGLE
! DUM =         DUMMY VARIABLE FOR THE CALL TO ARCINT
!
!
! COMPUTE WEIGHT FUNCTIONS C1, C2, AND C3
!
      c1 = b2*b3
      c2 = b3*b1
      c3 = b1*b2
      sum = c1 + c2 + c3
      if ( sum>0. ) then
!
! NORMALIZE C1, C2, AND C3
!
         c1 = c1/sum
         c2 = c2/sum
         c3 = c3/sum
!
! COMPUTE (U1,U2,U3) AND (U1N,U2N,U3N)
!
         u1n = 0.
         u2n = 0.
         u3n = 0.
         do i = 1 , 3
            u1(i) = (b2*v2(i)+b3*v3(i))/(b2+b3)
            u2(i) = (b3*v3(i)+b1*v1(i))/(b3+b1)
            u3(i) = (b1*v1(i)+b2*v2(i))/(b1+b2)
            u1n = u1n + u1(i)*u1(i)
            u2n = u2n + u2(i)*u2(i)
            u3n = u3n + u3(i)*u3(i)
         enddo
!
! COMPUTE Q1, Q2, AND Q3
!
         u1n = sqrt(u1n)
         u2n = sqrt(u2n)
         u3n = sqrt(u3n)
         do i = 1 , 3
            q1(i) = u1(i)/u1n
            q2(i) = u2(i)/u2n
            q3(i) = u3(i)/u3n
         enddo
!
! COMPUTE INTERPOLATED VALUE (VAL) AT P BY LOOPING ON
!   TRIANGLE SIDES
!
         val = 0.
!
! CONTRIBUTION FROM SIDE OPPOSITE V1 --
!
!   COMPUTE VALUE AND GRADIENT AT Q1 BY INTERPOLATING
!     BETWEEN V2 AND V3
!
         call arcint(q1,v2,v3,w2,w3,g2,g3,w,g,dum)
!
!   ADD IN THE CONTRIBUTION
!
         val = val + c1*(w+b1*b1*(3.-2.*b1)*(w1-w)+b1*(1.-b1)           &
     &         *(b1*(g1(1)*u1(1)+g1(2)*u1(2)+g1(3)*u1(3))+(1.-b1)       &
     &         *(g(1)*v1(1)+g(2)*v1(2)+g(3)*v1(3))/u1n))
!
! CONTRIBUTION FROM SIDE OPPOSITE V2 --
!
!   COMPUTE VALUE AND GRADIENT AT Q2 BY INTERPOLATING
!     BETWEEN V3 AND V1
!
         call arcint(q2,v3,v1,w3,w1,g3,g1,w,g,dum)
!
!   ADD IN THE CONTRIBUTION
!
         val = val + c2*(w+b2*b2*(3.-2.*b2)*(w2-w)+b2*(1.-b2)           &
     &         *(b2*(g2(1)*u2(1)+g2(2)*u2(2)+g2(3)*u2(3))+(1.-b2)       &
     &         *(g(1)*v2(1)+g(2)*v2(2)+g(3)*v2(3))/u2n))
!
! CONTRIBUTION FROM SIDE OPPOSITE V3 --
!
!   COMPUTE INTERPOLATED VALUE AND GRADIENT AT Q3
!     BY INTERPOLATING BETWEEN V1 AND V2
!
         call arcint(q3,v1,v2,w1,w2,g1,g2,w,g,dum)
!
!   ADD IN THE FINAL CONTRIBUTION
!
         val = val + c3*(w+b3*b3*(3.-2.*b3)*(w3-w)+b3*(1.-b3)           &
     &         *(b3*(g3(1)*u3(1)+g3(2)*u3(2)+g3(3)*u3(3))+(1.-b3)       &
     &         *(g(1)*v3(1)+g(2)*v3(2)+g(3)*v3(3))/u3n))
         pw = val
         if ( iflag/=1 ) return
         goto 99999
      endif
!
! P COINCIDES WITH A VERTEX
!
      pw = b1*w1 + b2*w2 + b3*w3
      if ( iflag/=1 ) return
      do i = 1 , 3
         pg(i) = b1*g1(i) + b2*g2(i) + b3*g3(i)
      enddo
      return
!
! COMPUTE PG
!*** THIS AIN'T YET BEEN CODED ***
!
99999 end subroutine wval
end module alg623
