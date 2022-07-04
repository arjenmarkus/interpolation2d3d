module alg624
contains
!*==ADNODE.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
!     ALGORITHM 624 COLLECTED ALGORITHMS FROM ACM.
!     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 4,
!     DEC., 1984, P. 453.
      subroutine adnode(kk,x,y,iadj,iend,ier)
      implicit none
!*--ADNODE7
!*** Start of declarations inserted by SPAG
      !integer index
!*** End of declarations inserted by SPAG
      integer kk , iadj(*) , iend(kk) , ier
      real x(kk) , y(kk)
      !logical swptst
      !external index
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS NODE KK TO A TRIANGULATION OF A SET
! OF POINTS IN THE PLANE PRODUCING A NEW TRIANGULATION.  A
! SEQUENCE OF EDGE SWAPS IS THEN APPLIED TO THE MESH,
! RESULTING IN AN OPTIMAL TRIANGULATION.  ADNODE IS PART
! OF AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES ROUTINES
! TO INITIALIZE THE DATA STRUCTURE, PLOT THE MESH, AND
! DELETE ARCS.
!
! INPUT PARAMETERS -   KK - INDEX OF THE NODE TO BE ADDED
!                           TO THE MESH.  KK .GE. 4.
!
!                     X,Y - VECTORS OF COORDINATES OF THE
!                           NODES IN THE MESH.  (X(I),Y(I))
!                           DEFINES NODE I FOR I = 1,..,KK.
!
!                    IADJ - SET OF ADJACENCY LISTS OF NODES
!                           1,..,KK-1.
!
!                    IEND - POINTERS TO THE ENDS OF
!                           ADJACENCY LISTS IN IADJ FOR
!                           EACH NODE IN THE MESH.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! KK, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.
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
!                                SHIFTD, INDEX, SWPTST,
!                                SWAP
!
!***********************************************************
!
      integer k , km1 , i1 , i2 , i3 , indkf , indkl , nabor1 , io1 ,   &
     &        io2 , in1 , indk1 , ind2f , ind21
      real xk , yk
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
! XK,YK =    X(K), Y(K)
!
      ier = 0
      k = kk
!
! INITIALIZATION
!
      km1 = k - 1
      xk = x(k)
      yk = y(k)
!
! ADD NODE K TO THE MESH
!
      call trfind(km1,xk,yk,x,y,iadj,iend,i1,i2,i3)
      if ( i1==0 ) then
!
! ALL NODES ARE COLLINEAR
!
         ier = 1
         goto 99999
      else
         if ( i3==0 ) call bdyadd(k,i1,i2,iadj,iend)
         if ( i3/=0 ) call intadd(k,i1,i2,i3,iadj,iend)
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
      if ( swptst(in1,k,io1,io2,x,y) ) then
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
!*==AREA.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      function area(x,y,nb,nodes)
      implicit none
!*--AREA161
!*** Start of declarations inserted by SPAG
      real area
!*** End of declarations inserted by SPAG
      integer nb , nodes(nb)
      real x(*) , y(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A SEQUENCE OF NB POINTS IN THE PLANE, THIS
! FUNCTION COMPUTES THE AREA BOUNDED BY THE CLOSED POLY-
! GONAL CURVE WHICH PASSES THROUGH THE POINTS IN THE
! SPECIFIED ORDER.  EACH SIMPLE CLOSED CURVE IS POSITIVELY
! ORIENTED (BOUNDS POSITIVE AREA) IF AND ONLY IF THE POINTS
! ARE SPECIFIED IN COUNTERCLOCKWISE ORDER.  THE LAST POINT
! OF THE CURVE IS TAKEN TO BE THE FIRST POINT SPECIFIED, AND
! THUS THIS POINT NEED NOT BE SPECIFIED TWICE.  HOWEVER, ANY
! POINT MAY BE SPECIFIED MORE THAN ONCE IN ORDER TO DEFINE A
! MULTIPLY CONNECTED DOMAIN.
!   THE AREA OF A TRIANGULATION MAY BE COMPUTED BY CALLING
! AREA WITH VALUES OF NB AND NODES DETERMINED BY SUBROUTINE
! BNODES.
!
! INPUT PARAMETERS -   X,Y - N-VECTORS OF COORDINATES OF
!                            POINTS IN THE PLANE FOR N .GE.
!                            NB.  NODE I HAS COORDINATES
!                            (X(I),Y(I)) FOR I = 1, 2, ...,
!                            N.
!
!                       NB - LENGTH OF NODES.
!
!                    NODES - VECTOR OF NODE INDICES IN THE
!                            RANGE 1 TO N DEFINING THE
!                            POLYGONAL CURVE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER -  AREA - SIGNED AREA BOUNDED BY THE
!                            POLYGONAL CURVE DEFINED
!                            ABOVE.
!
! MODULES REFERENCED BY AREA - NONE
!
!***********************************************************
!
      integer nnb , nd , i
      real a , x0 , y0 , dx1 , dy1 , dx2 , dy2
!
! LOCAL PARAMETERS -
!
! NNB =     LOCAL COPY OF NB
! ND =      ELEMENT OF NODES
! I =       DO-LOOP AND NODES INDEX
! A =       PARTIAL SUM OF SIGNED (AND DOUBLED) TRIANGLE
!             AREAS
! X0,Y0 =   X(NODES(1)), Y(NODES(1))
! DX1,DY1 = COMPONENTS OF THE VECTOR NODES(1)-NODES(I) FOR
!             I = 2, 3, ..., NB-1
! DX2,DY2 = COMPONENTS OF THE VECTOR NODES(1)-NODES(I) FOR
!             I = 3, 4, ..., NB
!
      nnb = nb
      a = 0.
      if ( nnb>=3 ) then
!
! INITIALIZATION
!
         nd = nodes(1)
         x0 = x(nd)
         y0 = y(nd)
         nd = nodes(2)
         dx1 = x(nd) - x0
         dy1 = y(nd) - y0
!
! LOOP ON TRIANGLES (NODES(1),NODES(I),NODES(I+1)),
!   I = 2, 3, ..., NB-1, ADDING TWICE THEIR SIGNED
!   AREAS TO A
!
         do i = 3 , nnb
            nd = nodes(i)
            dx2 = x(nd) - x0
            dy2 = y(nd) - y0
            a = a + dx1*dy2 - dx2*dy1
            dx1 = dx2
            dy1 = dy2
         enddo
      endif
!
! A CONTAINS TWICE THE SIGNED AREA OF THE REGION
!
      area = a/2.
      end function area
!*==BDYADD.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine bdyadd(kk,i1,i2,iadj,iend)
      implicit none
!*--BDYADD260
      integer kk , i1 , i2 , iadj(*) , iend(kk)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS A BOUNDARY NODE TO A TRIANGULATION
! OF A SET OF KK-1 POINTS IN THE PLANE.  IADJ AND IEND ARE
! UPDATED WITH THE INSERTION OF NODE KK.
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
!*==BNODES.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine bnodes(n,iadj,iend,nb,na,nt,nodes)
      implicit none
!*--BNODES401
      integer n , iadj(*) , iend(n) , nb , na , nt , nodes(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N POINTS IN THE PLANE, THIS
! ROUTINE RETURNS A VECTOR CONTAINING THE INDICES, IN
! COUNTERCLOCKWISE ORDER, OF THE NODES ON THE BOUNDARY OF
! THE CONVEX HULL OF THE SET OF POINTS.
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
      integer nst , indl , k , n0 , indf
!
! LOCAL PARAMETERS -
!
! NST =  FIRST ELEMENT OF NODES -- ARBITRARILY CHOSEN
! INDL = IADJ INDEX OF THE LAST NEIGHBOR OF NST
! K =    NODES INDEX
! N0 =   BOUNDARY NODE TO BE ADDED TO NODES
! INDF = IADJ INDEX OF THE FIRST NEIGHBOR OF N0
!
! SET NST TO THE FIRST BOUNDARY NODE ENCOUNTERED
!
      nst = 1
      do
         indl = iend(nst)
         if ( iadj(indl)==0 ) then
!
! INITIALIZATION
!
            nodes(1) = nst
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
                  nt = 2*n - nb - 2
                  na = nt + n - 1
                  goto 99999
               else
                  k = k + 1
                  nodes(k) = n0
               endif
            enddo
         else
            nst = nst + 1
         endif
      enddo
99999 end subroutine bnodes
!*==DELETE.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine delete(nn,nout1,nout2,iadj,iend,ier)
      implicit none
!*--DELETE492
!*** Start of declarations inserted by SPAG
      !integer index
!*** End of declarations inserted by SPAG
      integer nn , nout1 , nout2 , iadj(*) , iend(nn) , ier
      !external index
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE DELETES A BOUNDARY EDGE FROM A TRIANGU-
! LATION OF A SET OF POINTS IN THE PLANE.  IT MAY BE NEC-
! ESSARY TO FORCE CERTAIN EDGES TO BE PRESENT BEFORE CALL-
! ING DELETE (SEE SUBROUTINE EDGE).  NOTE THAT SUBROUTINES
! EDGE, TRFIND, AND THE ROUTINES WHICH CALL TRFIND (ADNODE,
! UNIF, INTRC1, AND INTRC0) SHOULD NOT BE CALLED FOLLOWING
! A DELETION.
!
! INPUT PARAMETERS -    NN - NUMBER OF NODES IN THE TRIAN-
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
      integer n , iout1 , iout2 , io1 , io2 , ind12 , ind21 , itemp ,   &
     &        ind1f , ind1l , ind2f , ind2l , newbd , indnf , indnl ,   &
     &        indn0 , indfp2 , indlm3 , nf , nl , i , imax
!
! LOCAL PARAMETERS -
!
! N =           LOCAL COPY OF NN
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
      n = nn
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
! ONE OF THE VERTICES HAS ONLY TWO NONZERO NEIGHBORS.  THE
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
               nl = iend(n)
               if ( nf<=nl ) call shiftd(nf,nl,-1,iadj)
               do i = newbd , n
                  iend(i) = iend(i) - 1
               enddo
               goto 50
            endif
!
! DELETE IO1 AS A NEIGHBOR OF IO2
!
            nf = ind21 + 1
            nl = iend(n)
            call shiftd(nf,nl,-1,iadj)
            do i = io2 , n
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
!*==EDGE.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine edge(in1,in2,x,y,lwk,iwk,iadj,iend,ier)
      implicit none
!*--EDGE780
!*** Start of declarations inserted by SPAG
      real xa , xb , xp , ya , yb , yp
!*** End of declarations inserted by SPAG
      !logical swptst
      integer in1 , in2 , lwk , iwk(2,*) , iadj(*) , iend(*) , ier
      real x(*) , y(*)
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
! INPUT PARAMETERS - IN1,IN2 - INDICES (OF X AND Y) IN THE
!                              RANGE 1,...,N DEFINING A PAIR
!                              OF NODES TO BE CONNECTED BY
!                              AN ARC.
!
!                        X,Y - N-VECTORS CONTAINING CARTE-
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
      real x1 , y1 , x2 , y2 , x0 , y0
      logical swp , left
!
! LOCAL PARAMETERS -
!
! N1,N2 =   LOCAL COPIES OF IN1 AND IN2 OR NODES OPPOSITE AN
!             ARC IO1-IO2 TO BE TESTED FOR A SWAP IN THE
!             OPTIMIZATION LOOPS
! IWEND =   INPUT OR OUTPUT VALUE OF LWK
! IWL =     IWK (COLUMN) INDEX OF THE LAST (RIGHTMOST) ARC
!             WHICH INTERSECTS IN1->IN2
! INDF =    IADJ INDEX OF THE FIRST NEIGHBOR OF IN1 OR IO1
! INDX =    IADJ INDEX OF A NEIGHBOR OF IN1, NL, OR IO1
! N1LST =   LAST NEIGHBOR OF IN1
! NL,NR =   ENDPOINTS OF AN ARC WHICH INTERSECTS IN1-IN2
!             WITH NL LEFT IN1->IN2
! NEXT =    NODE OPPOSITE NL->NR
! IWF =     IWK (COLUMN) INDEX OF THE FIRST (LEFTMOST) ARC
!             WHICH INTERSECTS IN1->IN2
! LFT =     FLAG USED TO DETERMINE IF A SWAP RESULTS IN THE
!             NEW ARC INTERSECTING IN1-IN2 -- LFT = 0 IFF
!             N0 = IN1, LFT = -1 IMPLIES N0 LEFT IN1->IN2,
!             AND LFT = 1 IMPLIES N0 LEFT IN2->IN1
! N0 =      NODE OPPOSITE NR->NL
! IWC =     IWK INDEX BETWEEN IWF AND IWL -- NL->NR IS
!             STORED IN IWK(1,IWC)->IWK(2,IWC)
! IWCP1 =   IWC + 1
! IWCM1 =   IWC - 1
! I =       DO-LOOP INDEX AND COLUMN INDEX FOR IWK
! IO1,IO2 = ENDPOINTS OF AN ARC TO BE TESTED FOR A SWAP IN
!             THE OPTIMIZATION LOOPS
! INDL =    IADJ INDEX OF THE LAST NEIGHBOR OF IO1
! X1,Y1 =   COORDINATES OF IN1
! X2,Y2 =   COORDINATES OF IN2
! X0,Y0 =   COORDINATES OF N0
! SWP =     FLAG SET TO .TRUE. IFF A SWAP OCCURS IN AN OPTI-
!             MIZATION LOOP
! LEFT =    STATEMENT FUNCTION WHICH RETURNS THE VALUE
!             .TRUE. IFF (XP,YP) IS ON OR TO THE LEFT OF THE
!             VECTOR (XA,YA)->(XB,YB)
!
      left(xa,ya,xb,yb,xp,yp) = (xb-xa)*(yp-ya)>=(xp-xa)*(yb-ya)
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
         x2 = x(n2)
         y2 = y(n2)
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
            do while ( .not.(left(x1,y1,x2,y2,x(nl),y(nl))) )
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
            if ( left(x2,y2,x1,y1,x(nr),y(nr)) ) then
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
         if ( left(x1,y1,x2,y2,x(next),y(next)) ) then
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
            if ( .not.left(x(nl),y(nl),x0,y0,x(next),y(next)) ) goto 400
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
               if ( .not.left(x0,y0,x(nr),y(nr),x(next),y(next)) )      &
     &              goto 400
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
            if ( left(x0,y0,x(nr),y(nr),x(next),y(next)) ) then
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
               elseif ( left(x(nl),y(nl),x0,y0,x(next),y(next)) ) then
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
               if ( swptst(n1,n2,io1,io2,x,y) ) then
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
         if ( .not.left(x(nl),y(nl),x0,y0,x2,y2) ) goto 200
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
         if ( left(x0,y0,x(nr),y(nr),x2,y2) ) then
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
            if ( swptst(n1,n2,io1,io2,x,y) ) then
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
!*==GETNP.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine getnp(x,y,iadj,iend,l,npts,ds,ier)
      implicit none
!*--GETNP1358
      integer iadj(*) , iend(*) , l , npts(*) , ier
      real x(*) , y(*) , ds
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF N NODES AND AN ARRAY
! NPTS CONTAINING THE INDICES OF L-1 NODES ORDERED BY
! EUCLIDEAN DISTANCE FROM NPTS(1), THIS SUBROUTINE SETS
! NPTS(L) TO THE INDEX OF THE NEXT NODE IN THE SEQUENCE --
! THE NODE, OTHER THAN NPTS(1),...,NPTS(L-1), WHICH IS
! CLOSEST TO NPTS(1).  THUS, THE ORDERED SEQUENCE OF K
! CLOSEST NODES TO N1 (INCLUDING N1) MAY BE DETERMINED BY
! K-1 CALLS TO GETNP WITH NPTS(1) = N1 AND L = 2,3,...,K
! FOR K .GE. 2.
!   THE ALGORITHM USES THE FACT THAT, IN A THIESSEN TRIAN-
! GULATION, THE K-TH CLOSEST NODE TO A GIVEN NODE N1 IS A
! NEIGHBOR OF ONE OF THE K-1 CLOSEST NODES TO N1.
!
! INPUT PARAMETERS - X,Y - VECTORS OF LENGTH N CONTAINING
!                          THE CARTESIAN COORDINATES OF THE
!                          NODES.
!
!                   IADJ - SET OF ADJACENCY LISTS OF NODES
!                          IN THE TRIANGULATION.
!
!                   IEND - POINTERS TO THE ENDS OF ADJACENCY
!                          LISTS FOR EACH NODE IN THE TRI-
!                          ANGULATION.
!
!                      L - NUMBER OF NODES IN THE SEQUENCE
!                          ON OUTPUT.  2 .LE. L .LE. N.
!
!                   NPTS - ARRAY OF LENGTH .GE. L CONTAIN-
!                          ING THE INDICES OF THE L-1 CLOS-
!                          EST NODES TO NPTS(1) IN THE FIRST
!                          L-1 LOCATIONS.
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
!                       DS - SQUARED EUCLIDEAN DISTANCE BE-
!                            TWEEN NPTS(1) AND NPTS(L)
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
      real x1 , y1 , dnp , dnb
!
! LOCAL PARAMETERS -
!
! LM1 =     L - 1
! N1 =      NPTS(1)
! I =       NPTS INDEX AND DO-LOOP INDEX
! NI =      NPTS(I)
! NP =      CANDIDATE FOR NPTS(L)
! INDF =    IADJ INDEX OF THE FIRST NEIGHBOR OF NI
! INDL =    IADJ INDEX OF THE LAST NEIGHBOR OF NI
! INDX =    IADJ INDEX IN THE RANGE INDF,...,INDL
! NB =      NEIGHBOR OF NI AND CANDIDATE FOR NP
! X1,Y1 =   COORDINATES OF N1
! DNP,DNB = SQUARED DISTANCES FROM N1 TO NP AND NB,
!             RESPECTIVELY
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
!
! MARK THE ELEMENTS OF NPTS
!
      do i = 1 , lm1
         ni = npts(i)
         iend(ni) = -iend(ni)
      enddo
!
! CANDIDATES FOR NP = NPTS(L) ARE THE UNMARKED NEIGHBORS
!   OF NODES IN NPTS.  NP=0 IS A FLAG TO SET NP TO THE
!   FIRST CANDIDATE ENCOUNTERED.
!
      np = 0
      dnp = 0.
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
            if ( nb/=0 ) then
               if ( iend(nb)>=0 ) then
!
! NB IS AN UNMARKED NEIGHBOR OF NI.  REPLACE NP IF NB IS
!   CLOSER TO N1 OR IS THE FIRST CANDIDATE ENCOUNTERED.
!
                  dnb = (x(nb)-x1)**2 + (y(nb)-y1)**2
                  if ( np==0 .or. dnb<dnp ) then
                     np = nb
                     dnp = dnb
                  endif
               endif
            endif
         enddo
      enddo
      npts(l) = np
      ds = dnp
!
! UNMARK THE ELEMENTS OF NPTS
!
      do i = 1 , lm1
         ni = npts(i)
         iend(ni) = -iend(ni)
      enddo
      return
99999 end subroutine getnp
!*==INDEX.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      integer function indexr(nvertx,nabor,iadj,iend)
      implicit none
!*--INDEX1508
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
!*==INTADD.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine intadd(kk,i1,i2,i3,iadj,iend)
      implicit none
!*--INTADD1567
      integer kk , i1 , i2 , i3 , iadj(*) , iend(kk)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE ADDS AN INTERIOR NODE TO A TRIANGULATION
! OF A SET OF KK-1 POINTS IN THE PLANE.  IADJ AND IEND ARE
! UPDATED WITH THE INSERTION OF NODE KK IN THE TRIANGLE
! WHOSE VERTICES ARE I1, I2, AND I3.
!
! INPUT PARAMETERS -        KK - INDEX OF NODE TO BE
!                                INSERTED.  KK .GE. 4.
!
!                     I1,I2,I3 - INDICES OF THE VERTICES OF
!                                A TRIANGLE CONTAINING NODE
!                                KK -- IN COUNTERCLOCKWISE
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
! ORDER THE VERTICES BY DECREASING MAGNITUDE.
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
!*==PERMUT.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine permut(nn,ip,a)
      implicit none
!*--PERMUT1732
      integer nn , ip(nn)
      real a(nn)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE APPLIES A SET OF PERMUTATIONS TO A VECTOR.
!
! INPUT PARAMETERS - NN - LENGTH OF A AND IP.
!
!                    IP - VECTOR CONTAINING THE SEQUENCE OF
!                         INTEGERS 1,...,NN PERMUTED IN THE
!                         SAME FASHION THAT A IS TO BE PER-
!                         MUTED.
!
!                     A - VECTOR TO BE PERMUTED.
!
! NN AND IP ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - A - REORDERED VECTOR REFLECTING THE
!                         PERMUTATIONS DEFINED BY IP.
!
! MODULES REFERENCED BY PERMUT - NONE
!
!***********************************************************
!
      integer n , k , j , ipj
      real temp
!
! LOCAL PARAMETERS -
!
! N =    LOCAL COPY OF NN
! K =    INDEX FOR IP AND FOR THE FIRST ELEMENT OF A IN A
!          PERMUTATION
! J =    INDEX FOR IP AND A, J .GE. K
! IPJ =  IP(J)
! TEMP = TEMPORARY STORAGE FOR A(K)
!
      n = nn
      if ( n<2 ) return
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
               if ( k>n ) then
!
! ALL PERMUTATIONS HAVE BEEN APPLIED.  UNMARK IP.
!
                  do k = 1 , n
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
!*==QSORT.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine qsort(n,x,ind)
      implicit none
!*--QSORT1817
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
!*==REORDR.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine reordr(n,iflag,a,b,c,ind)
      implicit none
!*--REORDR2042
      integer n , iflag , ind(n)
      real a(n) , b(n) , c(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO
! REORDER THE REAL ARRAY A INTO INCREASING ORDER.  A RECORD
! OF THE PERMUTATIONS APPLIED TO A IS STORED IN IND, AND
! THESE PERMUTATIONS MAY BE APPLIED TO ONE OR TWO ADDITIONAL
! VECTORS BY THIS ROUTINE.  ANY OTHER VECTOR V MAY BE PER-
! MUTED IN THE SAME FASHION BY CALLING SUBROUTINE PERMUT
! WITH N, IND, AND V AS PARAMETERS.
!   A SET OF NODES (X(I),Y(I)) AND DATA VALUES Z(I) MAY BE
! PREPROCESSED BY REORDR FOR INCREASED EFFICIENCY IN THE
! TRIANGULATION ROUTINE TRMESH.  EFFICIENCY IS INCREASED BY
! A FACTOR OF APPROXIMATELY SQRT(N)/6 FOR RANDOMLY DISTRIB-
! UTED NODES, AND THE PREPROCESSING IS ALSO USEFUL FOR
! DETECTING DUPLICATE NODES.  EITHER X OR Y MAY BE USED AS
! THE SORT KEY (ASSOCIATED WITH A).
!
! INPUT PARAMETERS - N - NUMBER OF NODES.
!
!                IFLAG - NUMBER OF VECTORS TO BE PERMUTED.
!                        IFLAG .LE. 0 IF A, B, AND C ARE TO
!                                     REMAIN UNALTERED.
!                        IFLAG .EQ. 1 IF ONLY A IS TO BE
!                                     PERMUTED.
!                        IFLAG .EQ. 2 IF A AND B ARE TO BE
!                                     PERMUTED.
!                        IFLAG .GE. 3 IF A, B, AND C ARE TO
!                                     BE PERMUTED.
!
!                A,B,C - VECTORS OF LENGTH N TO BE SORTED
!                        (ON THE COMPONENTS OF A), OR DUMMY
!                        PARAMETERS, DEPENDING ON IFLAG.
!
!                  IND - VECTOR OF LENGTH .GE. N.
!
! N, IFLAG, AND ANY DUMMY PARAMETERS ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - A,B,C - SORTED OR UNALTERED VECTORS.
!
!                       IND - SEQUENCE OF INDICES 1,...,N
!                             PERMUTED IN THE SAME FASHION
!                             AS THE REAL VECTORS.  THUS,
!                             THE ORDERING MAY BE APPLIED TO
!                             A REAL VECTOR V AND STORED IN
!                             W BY SETTING W(I) = V(IND(I)),
!                             OR V MAY BE OVERWRITTEN WITH
!                             THE ORDERING BY A CALL TO PER-
!                             MUT.
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
      end subroutine reordr
!*==SHIFTD.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine shiftd(nfrst,nlast,kk,iarr)
      implicit none
!*--SHIFTD2124
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
!*==SWAP.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine swap(nin1,nin2,nout1,nout2,iadj,iend)
      implicit none
!*--SWAP2244
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
!*==SWPTST.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      logical function swptst(in1,in2,io1,io2,x,y)
      implicit none
!*--SWPTST2447
      integer in1 , in2 , io1 , io2
      real x(*) , y(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION DECIDES WHETHER OR NOT TO REPLACE A
! DIAGONAL ARC IN A QUADRILATERAL WITH THE OTHER DIAGONAL.
! THE DETERMINATION IS BASED ON THE SIZES OF THE ANGLES
! CONTAINED IN THE 2 TRIANGLES DEFINED BY THE DIAGONAL.
! THE DIAGONAL IS CHOSEN TO MAXIMIZE THE SMALLEST OF THE
! SIX ANGLES OVER THE TWO PAIRS OF TRIANGLES.
!
! INPUT PARAMETERS -  IN1,IN2,IO1,IO2 - NODE INDICES OF THE
!                              FOUR POINTS DEFINING THE
!                              QUADRILATERAL.  IO1 AND IO2
!                              ARE CURRENTLY CONNECTED BY A
!                              DIAGONAL ARC.  THIS ARC
!                              SHOULD BE REPLACED BY AN ARC
!                              CONNECTING IN1, IN2 IF THE
!                              DECISION IS MADE TO SWAP.
!                              IN1,IO1,IO2 MUST BE IN
!                              COUNTERCLOCKWISE ORDER.
!
!                        X,Y - VECTORS OF NODAL COORDINATES.
!                              (X(I),Y(I)) ARE THE COORD-
!                              INATES OF NODE I FOR I = IN1,
!                              IN2, IO1, OR IO2.
!
! NONE OF THE INPUT PARAMETERS ARE ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER -  SWPTST - .TRUE. IFF THE ARC CONNECTING
!                              IO1 AND IO2 IS TO BE REPLACED
!
! MODULES REFERENCED BY SWPTST - NONE
!
!***********************************************************
!
      real dx11 , dx12 , dx22 , dx21 , dy11 , dy12 , dy22 , dy21 ,      &
     &     sin1 , sin2 , cos1 , cos2 , sin12
!
! LOCAL PARAMETERS -
!
! DX11,DY11 = X,Y COORDINATES OF THE VECTOR IN1-IO1
! DX12,DY12 = X,Y COORDINATES OF THE VECTOR IN1-IO2
! DX22,DY22 = X,Y COORDINATES OF THE VECTOR IN2-IO2
! DX21,DY21 = X,Y COORDINATES OF THE VECTOR IN2-IO1
! SIN1 =      CROSS PRODUCT OF THE VECTORS IN1-IO1 AND
!               IN1-IO2 -- PROPORTIONAL TO SIN(T1) WHERE T1
!               IS THE ANGLE AT IN1 FORMED BY THE VECTORS
! COS1 =      INNER PRODUCT OF THE VECTORS IN1-IO1 AND
!               IN1-IO2 -- PROPORTIONAL TO COS(T1)
! SIN2 =      CROSS PRODUCT OF THE VECTORS IN2-IO2 AND
!               IN2-IO1 -- PROPORTIONAL TO SIN(T2) WHERE T2
!               IS THE ANGLE AT IN2 FORMED BY THE VECTORS
! COS2 =      INNER PRODUCT OF THE VECTORS IN2-IO2 AND
!               IN2-IO1 -- PROPORTIONAL TO COS(T2)
! SIN12 =     SIN1*COS2 + COS1*SIN2 -- PROPORTIONAL TO
!               SIN(T1+T2)
!
      swptst = .false.
!
! COMPUTE THE VECTORS CONTAINING THE ANGLES T1, T2
!
      dx11 = x(io1) - x(in1)
      dx12 = x(io2) - x(in1)
      dx22 = x(io2) - x(in2)
      dx21 = x(io1) - x(in2)
!
      dy11 = y(io1) - y(in1)
      dy12 = y(io2) - y(in1)
      dy22 = y(io2) - y(in2)
      dy21 = y(io1) - y(in2)
!
! COMPUTE INNER PRODUCTS
!
      cos1 = dx11*dx12 + dy11*dy12
      cos2 = dx22*dx21 + dy22*dy21
!
! THE DIAGONALS SHOULD BE SWAPPED IFF (T1+T2) .GT. 180
!   DEGREES.  THE FOLLOWING TWO TESTS INSURE NUMERICAL
!   STABILITY.
!
      if ( cos1>=0. .and. cos2>=0. ) return
      if ( cos1>=0. .or. cos2>=0. ) then
!
! COMPUTE VECTOR CROSS PRODUCTS
!
         sin1 = dx11*dy12 - dx12*dy11
         sin2 = dx22*dy21 - dx21*dy22
         sin12 = sin1*cos2 + cos1*sin2
         if ( sin12>=0. ) return
      endif
      swptst = .true.
      end function swptst
!*==TRFIND.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine trfind(nst,px,py,x,y,iadj,iend,i1,i2,i3)
      implicit none
!*--TRFIND2549
!*** Start of declarations inserted by SPAG
      real x0 , x1 , x2 , y0 , y1 , y2
!*** End of declarations inserted by SPAG
      integer nst , iadj(*) , iend(*) , i1 , i2 , i3
      real px , py , x(*) , y(*)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE LOCATES A POINT P IN A THIESSEN TRIANGU-
! LATION, RETURNING THE VERTEX INDICES OF A TRIANGLE WHICH
! CONTAINS P.  TRFIND IS PART OF AN INTERPOLATION PACKAGE
! WHICH PROVIDES SUBROUTINES FOR CREATING THE MESH.
!
! INPUT PARAMETERS -    NST - INDEX OF NODE AT WHICH TRFIND
!                             BEGINS SEARCH.  SEARCH TIME
!                             DEPENDS ON THE PROXIMITY OF
!                             NST TO P.
!
!                     PX,PY - X AND Y-COORDINATES OF THE
!                             POINT TO BE LOCATED.
!
!                       X,Y - VECTORS OF COORDINATES OF
!                             NODES IN THE MESH.  (X(I),Y(I))
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
! OUTPUT PARAMETERS - I1,I2,I3 - VERTEX INDICES IN COUNTER-
!                                CLOCKWISE ORDER - VERTICES
!                                OF A TRIANGLE CONTAINING P
!                                IF P IS AN INTERIOR NODE.
!                                IF P IS OUTSIDE OF THE
!                                BOUNDARY OF THE MESH, I1
!                                AND I2 ARE THE FIRST (RIGHT
!                                -MOST) AND LAST (LEFTMOST)
!                                NODES WHICH ARE VISIBLE
!                                FROM P, AND I3 = 0.  IF P
!                                AND ALL OF THE NODES LIE ON
!                                A SINGLE LINE THEN I1 = I2
!                                = I3 = 0.
!
! MODULES REFERENCED BY TRFIND - NONE
!
! INTRINSIC FUNCTION CALLED BY TRFIND - MAX0
!
!***********************************************************
!
      integer n0 , n1 , n2 , n3 , n4 , indx , ind , nf , nl , next
      real xp , yp
      logical left
!
! LOCAL PARAMETERS -
!
! XP,YP =     LOCAL VARIABLES CONTAINING PX AND PY
! N0,N1,N2 =  NODES IN COUNTERCLOCKWISE ORDER DEFINING A
!               CONE (WITH VERTEX N0) CONTAINING P
! N3,N4 =     NODES OPPOSITE N1-N2 AND N2-N1, RESPECTIVELY
! INDX,IND =  INDICES FOR IADJ
! NF,NL =     FIRST AND LAST NEIGHBORS OF N0 IN IADJ, OR
!               FIRST (RIGHTMOST) AND LAST (LEFTMOST) NODES
!               VISIBLE FROM P WHEN P IS OUTSIDE THE
!               BOUNDARY
! NEXT =      CANDIDATE FOR I1 OR I2 WHEN P IS OUTSIDE OF
!               THE BOUNDARY
! LEFT =      STATEMENT FUNCTION WHICH COMPUTES THE SIGN OF
!               A CROSS PRODUCT (Z-COMPONENT).  LEFT(X1,...,
!               Y0) = .TRUE. IFF (X0,Y0) IS ON OR TO THE
!               LEFT OF THE VECTOR FROM (X1,Y1) TO (X2,Y2).
!
      left(x1,y1,x2,y2,x0,y0) = (x2-x1)*(y0-y1)>=(x0-x1)*(y2-y1)
      xp = px
      yp = py
!
! INITIALIZE VARIABLES AND FIND A CONE CONTAINING P
!
      n0 = max0(nst,1)
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
         do while ( .not.(left(x(n0),y(n0),x(n1),y(n1),xp,yp)) )
            indx = indx + 1
            n1 = iadj(indx)
            if ( n1==nl ) then
!
! P IS BETWEEN ARCS N0-N1 AND N0-NF
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
         if ( .not.(left(x(n0),y(n0),x(nf),y(nf),xp,yp)) ) then
!
! P IS OUTSIDE THE BOUNDARY
!
            nl = n0
            goto 300
         elseif ( .not.(left(x(nl),y(nl),x(n0),y(n0),xp,yp)) ) then
!
! P IS OUTSIDE THE BOUNDARY AND N0 IS THE RIGHTMOST
!   VISIBLE BOUNDARY NODE
!
            i1 = n0
            goto 400
         endif
      endif
      do
!
! P IS TO THE LEFT OF ARC N0-N1.  INITIALIZE N2 TO THE NEXT
!   NEIGHBOR OF N0.
!
         indx = indx + 1
         n2 = iadj(indx)
         if ( .not.left(x(n0),y(n0),x(n2),y(n2),xp,yp) ) then
            n3 = n0
            exit
         else
            n1 = n2
            if ( n1==nl ) then
               if ( .not.left(x(n0),y(n0),x(nf),y(nf),xp,yp) ) then
                  n2 = nf
                  n3 = n0
                  exit
               elseif ( xp==x(n0) .and. yp==y(n0) ) then
!
! P IS TO THE RIGHT OF N1-N0, OR P=N0.  SET N0 TO N1 AND
!   START OVER.
!
                  n0 = n1
                  goto 100
               else
!
! P IS LEFT OF OR ON ARCS N0-NB FOR ALL NEIGHBORS NB
!   OF N0.
! ALL POINTS ARE COLLINEAR IFF P IS LEFT OF NB-N0 FOR
!   ALL NEIGHBORS NB OF N0.  SEARCH THE NEIGHBORS OF N0
!   IN REVERSE ORDER.  NOTE -- N1 = NL AND INDX POINTS TO
!   NL.
!
                  do while ( left(x(n1),y(n1),x(n0),y(n0),xp,yp) )
                     if ( n1==nf ) then
!
! ALL POINTS ARE COLLINEAR
!
                        i1 = 0
                        i2 = 0
                        i3 = 0
                        goto 99999
                     else
                        indx = indx - 1
                        n1 = iadj(indx)
                     endif
                  enddo
                  n0 = n1
                  goto 100
               endif
            endif
         endif
      enddo
 200  do while ( .not.(left(x(n1),y(n1),x(n2),y(n2),xp,yp)) )
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
!
! P IS OUTSIDE THE BOUNDARY
!
                  nf = n2
                  nl = n1
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
! DEFINE A NEW ARC N1-N2 WHICH INTERSECTS THE LINE
!   SEGMENT N0-P
!
         if ( left(x(n0),y(n0),x(n4),y(n4),xp,yp) ) then
            n3 = n1
            n1 = n4
         else
            n3 = n2
            n2 = n4
         endif
      enddo
!
! P IS IN THE TRIANGLE (N1,N2,N3) AND NOT ON N2-N3.  IF
!   N3-N1 OR N1-N2 IS A BOUNDARY ARC CONTAINING P, TREAT P
!   AS EXTERIOR.
!
      indx = iend(n1)
      if ( iadj(indx)==0 ) then
!
! N1 IS A BOUNDARY NODE.  N3-N1 IS A BOUNDARY ARC IFF N3
!   IS THE LAST NONZERO NEIGHBOR OF N1.
!
         if ( n3==iadj(indx-1) ) then
!
! N3-N1 IS A BOUNDARY ARC
!
            if ( left(x(n1),y(n1),x(n3),y(n3),xp,yp) ) then
!
! P LIES ON N1-N3
!
               i1 = n1
               i2 = n3
               i3 = 0
               return
            endif
         endif
!
! N3-N1 IS NOT A BOUNDARY ARC CONTAINING P.  N1-N2 IS A
!   BOUNDARY ARC IFF N2 IS THE FIRST NEIGHBOR OF N1.
!
         indx = 1
         if ( n1/=1 ) indx = iend(n1-1) + 1
         if ( n2==iadj(indx) ) then
!
! N1-N2 IS A BOUNDARY ARC
!
            if ( left(x(n2),y(n2),x(n1),y(n1),xp,yp) ) then
!
! P LIES ON N1-N2
!
               i1 = n2
               i2 = n1
               i3 = 0
               return
            endif
         endif
      endif
!
! P DOES NOT LIE ON A BOUNDARY ARC.
!
      i1 = n1
      i2 = n2
      i3 = n3
      return
 300  do
!
! NF AND NL ARE ADJACENT BOUNDARY NODES WHICH ARE VISIBLE
!   FROM P.  FIND THE FIRST VISIBLE BOUNDARY NODE.
! SET NEXT TO THE FIRST NEIGHBOR OF NF.
!
         indx = 1
         if ( nf/=1 ) indx = iend(nf-1) + 1
         next = iadj(indx)
         if ( left(x(nf),y(nf),x(next),y(next),xp,yp) ) then
!
! NF IS THE FIRST (RIGHTMOST) VISIBLE BOUNDARY NODE
!
            i1 = nf
            exit
         else
            nf = next
         endif
      enddo
 400  do
!
! FIND THE LAST VISIBLE BOUNDARY NODE.  NL IS THE FIRST
!   CANDIDATE FOR I2.
! SET NEXT TO THE LAST NEIGHBOR OF NL.
!
         indx = iend(nl) - 1
         next = iadj(indx)
         if ( left(x(next),y(next),x(nl),y(nl),xp,yp) ) exit
         nl = next
      enddo
!
! NL IS THE LAST (LEFTMOST) VISIBLE BOUNDARY NODE
!
      i2 = nl
      i3 = 0
      return
99999 end subroutine trfind
!*==TRMESH.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine trmesh(n,x,y,iadj,iend,ier)
      implicit none
!*--TRMESH2877
      integer n , iadj(*) , iend(n) , ier
      real x(n) , y(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE CREATES A THIESSEN TRIANGULATION OF N
! ARBITRARILY SPACED POINTS IN THE PLANE REFERRED TO AS
! NODES.  THE TRIANGULATION IS OPTIMAL IN THE SENSE THAT IT
! IS AS NEARLY EQUIANGULAR AS POSSIBLE.  TRMESH IS PART OF
! AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES SUBROUTINES
! TO REORDER THE NODES, ADD A NEW NODE, DELETE AN ARC, PLOT
! THE MESH, AND PRINT THE DATA STRUCTURE.
!   UNLESS THE NODES ARE ALREADY ORDERED IN SOME REASONABLE
! FASHION, THEY SHOULD BE REORDERED BY SUBROUTINE REORDR FOR
! INCREASED EFFICIENCY BEFORE CALLING TRMESH.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            N .GE. 3.
!
!                      X,Y - N-VECTORS OF COORDINATES.
!                            (X(I),Y(I)) DEFINES NODE I.
!
!                     IADJ - VECTOR OF LENGTH .GE. 6*N-9.
!
!                     IEND - VECTOR OF LENGTH .GE. N.
!
! N, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - IADJ - ADJACENCY LISTS OF NEIGHBORS IN
!                            COUNTERCLOCKWISE ORDER.  THE
!                            LIST FOR NODE I+1 FOLLOWS THAT
!                            FOR NODE I WHERE X AND Y DEFINE
!                            THE ORDER.  THE VALUE 0 DENOTES
!                            THE BOUNDARY (OR A PSEUDO-NODE
!                            AT INFINITY) AND IS ALWAYS THE
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
!                            IEND IS UNCHANGED IF IER = 1.
!                            IF IER = 2 IEND CONTAINS THE
!                            INDICES OF A SEQUENCE OF N
!                            NODES ORDERED FROM LEFT TO
!                            RIGHT WHERE LEFT AND RIGHT ARE
!                            DEFINED BY ASSUMING NODE 1 IS
!                            TO THE LEFT OF NODE 2.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE
!                                    ENCOUNTERED.
!                            IER = 1 IF N .LT. 3.
!                            IER = 2 IF N .GE. 3 AND ALL
!                                    NODES ARE COLLINEAR.
!
! MODULES REFERENCED BY TRMESH - SHIFTD, ADNODE, TRFIND,
!                                INTADD, BDYADD, SWPTST,
!                                SWAP, INDEX
!
!***********************************************************
!
      integer nn , k , km1 , nl , nr , ind , indx , n0 , itemp , ierr , &
     &        km1d2 , kmi , i , kmin
      real xl , yl , xr , yr , dxr , dyr , xk , yk , dxk , dyk , cprod ,&
     &     sprod
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! K =           NODE (INDEX) TO BE INSERTED INTO IEND
! KM1 =         K-1 - (VARIABLE) LENGTH OF IEND
! NL,NR =       IEND(1), IEND(KM1) -- LEFTMOST AND RIGHTMOST
!                 NODES IN IEND AS VIEWED FROM THE RIGHT OF
!                 1-2 WHEN IEND CONTAINS THE INITIAL ORDERED
!                 SET OF NODAL INDICES
! XL,YL,XR,YR = X AND Y COORDINATES OF NL AND NR
! DXR,DYR =     XR-XL, YR-YL
! XK,YK =       X AND Y COORDINATES OF NODE K
! DXK,DYK =     XK-XL, YK-YL
! CPROD =       VECTOR CROSS PRODUCT OF NL-NR AND NL-K --
!                 USED TO DETERMINE THE POSITION OF NODE K
!                 WITH RESPECT TO THE LINE DEFINED BY THE
!                 NODES IN IEND
! SPROD =       SCALAR PRODUCT USED TO DETERMINE THE
!                 INTERVAL CONTAINING NODE K WHEN K IS ON
!                 THE LINE DEFINED BY THE NODES IN IEND
! IND,INDX =    INDICES FOR IEND AND IADJ, RESPECTIVELY
! N0,ITEMP =    TEMPORARY NODES (INDICES)
! IERR =        DUMMY PARAMETER FOR CALL TO ADNODE
! KM1D2,KMI,I = KM1/2, K-I, DO-LOOP INDEX -- USED IN IEND
!                 REORDERING LOOP
! KMIN =        FIRST NODE INDEX SENT TO ADNODE
!
      nn = n
      ier = 1
      if ( nn<3 ) return
      ier = 0
!
! INITIALIZE IEND, NL, NR, AND K
!
      iend(1) = 1
      iend(2) = 2
      xl = x(1)
      yl = y(1)
      xr = x(2)
      yr = y(2)
      k = 2
!
! BEGIN LOOP ON NODES 3,4,...
!
 100  dxr = xr - xl
      dyr = yr - yl
!
! NEXT LOOP BEGINS HERE IF NL AND NR ARE UNCHANGED
!
      do while ( k/=nn )
         km1 = k
         k = km1 + 1
         xk = x(k)
         yk = y(k)
         dxk = xk - xl
         dyk = yk - yl
         cprod = dxr*dyk - dxk*dyr
         if ( cprod>0. ) then
!
! NODE K IS TO THE LEFT OF NL-NR.  REORDER IEND SO THAT NL
!   IS THE LEFTMOST NODE AS VIEWED FROM K.
!
            km1d2 = km1/2
            do i = 1 , km1d2
               kmi = k - i
               itemp = iend(i)
               iend(i) = iend(kmi)
               iend(kmi) = itemp
            enddo
            goto 200
         else
            if ( cprod<0. ) goto 200
!
! NODE K LIES ON THE LINE CONTAINING NODES 1,2,...,K-1.
!   SET SPROD TO (NL-NR,NL-K).
!
            sprod = dxr*dxk + dyr*dyk
            if ( sprod>0. ) then
!
! NODE K IS TO THE RIGHT OF NL.  FIND THE LEFTMOST NODE
!   N0 WHICH LIES TO THE RIGHT OF K.
!   SET SPROD TO (N0-NL,N0-K).
!
               do ind = 2 , km1
                  n0 = iend(ind)
                  sprod = (xl-x(n0))*(xk-x(n0)) + (yl-y(n0))*(yk-y(n0))
                  if ( sprod>=0. ) goto 120
               enddo
!
! NODE K IS TO THE RIGHT OF NR.  INSERT K AS THE LAST
!   (RIGHTMOST) NODE IN IEND AND SET NR TO K.
!
               iend(k) = k
               xr = xk
               yr = yk
            else
!
! NODE K IS TO THE LEFT OF NL.  INSERT K AS THE FIRST
!   (LEFTMOST) NODE IN IEND AND SET NL TO K.
!
               call shiftd(1,km1,1,iend)
               iend(1) = k
               xl = xk
               yl = yk
            endif
            goto 100
!
! NODE K LIES BETWEEN IEND(IND-1) AND IEND(IND).  INSERT K
!   IN IEND.
!
 120        call shiftd(ind,km1,1,iend)
            iend(ind) = k
         endif
      enddo
!
! ALL NODES ARE COLLINEAR
!
      ier = 2
      goto 99999
!
! NODE K IS TO THE RIGHT OF NL-NR.  CREATE A TRIANGULATION
!   CONSISTING OF NODES 1,2,...,K.
!
 200  nl = iend(1)
      nr = iend(km1)
!
! CREATE THE ADJACENCY LISTS FOR THE FIRST K-1 NODES.
!   INSERT NEIGHBORS IN REVERSE ORDER.  EACH NODE HAS FOUR
!   NEIGHBORS EXCEPT NL AND NR WHICH HAVE THREE.
!
      do ind = 1 , km1
         n0 = iend(ind)
         indx = 4*n0
         if ( n0>=nl ) indx = indx - 1
         if ( n0>=nr ) indx = indx - 1
         iadj(indx) = 0
         indx = indx - 1
         if ( ind<km1 ) iadj(indx) = iend(ind+1)
         if ( ind<km1 ) indx = indx - 1
         iadj(indx) = k
         if ( ind/=1 ) iadj(indx-1) = iend(ind-1)
      enddo
!
! CREATE THE ADJACENCY LIST FOR NODE K
!
      indx = 5*km1 - 1
      iadj(indx) = 0
      do ind = 1 , km1
         indx = indx - 1
         iadj(indx) = iend(ind)
      enddo
!
! REPLACE IEND ELEMENTS WITH POINTERS TO IADJ
!
      indx = 0
      do ind = 1 , km1
         indx = indx + 4
         if ( ind==nl .or. ind==nr ) indx = indx - 1
         iend(ind) = indx
      enddo
      indx = indx + k
      iend(k) = indx
!
! ADD THE REMAINING NODES TO THE TRIANGULATION
!
      if ( k==nn ) return
      kmin = k + 1
      do k = kmin , nn
         call adnode(k,x,y,iadj,iend,ierr)
      enddo
      return
99999 end subroutine trmesh
!*==TRPLOT.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine trplot(n,x,y,iadj,iend,ititle,nc,numbr,ier)
      implicit none
!*--TRPLOT3133
!*** Start of declarations inserted by SPAG
      integer nc
!*** End of declarations inserted by SPAG
      integer n , iadj(*) , iend(n) , ititle(*) , numbr , ier
      real x(n) , y(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS SUBROUTINE PLOTS THE ARCS OF A TRIANGULATION OF N
! NODES IN THE PLANE.  CARDS WITH C* IN THE FIRST TWO COL-
! UMNS MUST BE REPLACED WITH CALLS TO USER-SUPPLIED GRAPHICS
! SUBROUTINES IN ORDER TO GET ANY USE OUT OF THIS ROUTINE.
!
! INPUT PARAMETERS - N - NUMBER OF NODES IN THE TRIANGULA-
!                        TION.  N .GE. 3.
!
!                  X,Y - CARTESIAN COORDINATES OF THE NODES.
!
!            IADJ,IEND - TRIANGULATION DATA STRUCTURE (SEE
!                        SUBROUTINE TRMESH).
!
!               ITITLE - INTEGER ARRAY CONTAINING A LINE OF
!                        TEXT TO BE CENTERED ABOVE THE PLOT
!                        IF NC .GT. 0.  ITITLE MUST BE INIT-
!                        IALIZED WITH HOLLERITH CONSTANTS OR
!                        READ WITH AN A-FORMAT.  ITS DIMEN-
!                        SION DEPENDS ON NC AND THE NUMBER
!                        OF CHARACTERS STORED IN A COMPUTER
!                        WORD.
!
!                   NC - NUMBER OF CHARACTERS IN ITITLE.
!                        0 .LE. NC .LE. 40.  NO TITLE IS
!                        PLOTTED IF NC = 0.
!
!                NUMBR - OPTION INDICATOR.  IF NUMBR .NE. 0,
!                        THE NODAL INDICES ARE PLOTTTED NEXT
!                        TO THE NODES.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - IER - ERROR INDICATOR
!                          IER = 0 IF NO ERRORS WERE ENCOUN-
!                                  TERED.
!                          IER = 1 IF N OR NC IS OUT OF
!                                  RANGE.
!                          IER = 2 IF THE NODES LIE ON A
!                                  HORIZONTAL OR VERTICAL
!                                  LINE.
!
! MODULES REFERENCED BY TRPLOT - NONE
!
! INTRINSIC FUNCTIONS CALLED BY TRPLOT - AMIN1, AMAX1, IABS
!
!***********************************************************
!
      integer nn , i , n0 , indf , indl , indx , n1 , new
      real xmin , xmax , ymin , ymax , dx , dy , x0 , y0
      logical st0
!
! LOCAL PARAMETERS -
!
! NN =        LOCAL COPY OF N
! I =         DO-LOOP INDEX
! N0 =        NODE WHICH IS TO BE CONNECTED TO ITS NEIGHBORS
!               WITH LINE SEGMENTS
! INDF =      IADJ INDEX OF THE FIRST NEIGHBOR OF N0
! INDL =      IADJ INDEX OF THE LAST NEIGHBOR OF N0
! INDX =      IADJ INDEX IN THE RANGE INDF,...,INDL
! N1 =        NEIGHBOR OF N0
! NEW =       NEIGHBOR OF N0 AND CANDIDATE FOR NEXT VALUE
!               OF N0
! XMIN,XMAX = MINIMUM AND MAXIMUM X COORDINATES
! YMIN,YMAX = MINIMUM AND MAXIMUM Y COORDINATES
! DX,DY =     XMAX-XMIN AND YMAX-YMIN (DATA SPACE DIMEN-
!               SIONS)
! X0,Y0 =     COORDINATES OF N0
! ST0 =       SWITCH USED TO ALTERNATE DIRECTION OF PEN
!               MOVEMENT -- TRUE IFF PEN STARTS AT N0
!
      nn = n
!
! CHECK FOR INVALID PARAMETERS
!
      if ( nn<3 .or. nc<0 .or. nc>40 ) then
!
! N OR NC IS OUT OF RANGE
!
         ier = 1
         return
      else
         ier = 0
!
! COMPUTE DATA SPACE DIMENSIONS
!
         xmin = x(1)
         ymin = y(1)
         xmax = xmin
         ymax = ymin
         do i = 1 , nn
            xmin = amin1(x(i),xmin)
            ymin = amin1(y(i),ymin)
            xmax = amax1(x(i),xmax)
            ymax = amax1(y(i),ymax)
         enddo
         dx = xmax - xmin
         dy = ymax - ymin
         if ( dx==0. .or. dy==0. ) then
!
! NODES ARE COLLINEAR
!
            ier = 2
            goto 99999
         else
!*  COMMANDS WHICH PERFORM THE FOLLOWING FUNCTIONS SHOULD
!*  BE INSERTED HERE --
!*    INITIALIZE THE PLOTTING ENVIRONMENT IF NECESSARY,
!*    COMPUTE PLOTTER SPACE DIMENSIONS,
!*    ESTABLISH A LINEAR MAPPING FROM THE DATA SPACE TO
!*      THE PLOTTER SPACE,
!*    OPTIONALLY, DRAW AND LABEL AXES, AND
!*    PLOT THE TITLE IF NC .NE. 0.
!
! INITIALIZE FOR LOOP ON NODES.  EACH NODE N0 IS TO BE CON-
!   NECTED TO ITS NEIGHBORS WHICH HAVE LARGER INDICES.  N0
!   IS THEN MARKED BY MAKING THE CORRESPONDING IEND ENTRY
!   NEGATIVE, AND THE SEARCH FOR THE NEXT UNMARKED NODE BE-
!   GINS WITH THE NEIGHBORS OF N0.
!
            n0 = 1
            indf = 1
         endif
      endif
!
! TOP OF LOOP -- SET INDL, X0, AND Y0, AND INITIALIZE ST0
!   AND INDX
!
 100  indl = iend(n0)
      x0 = x(n0)
      y0 = y(n0)
      st0 = .true.
      indx = indl
      if ( iadj(indx)==0 ) indx = indx - 1
!
!   LOOP ON NEIGHBORS OF N0 IN REVERSE ORDER
!
 200  n1 = iadj(indx)
!
!   CONNECT N0 AND N1 -- THE DIRECTION OF PEN MOVEMENT
!     ALTERNATES BETWEEN AWAY FROM N0 AND TOWARD N0 FOR
!     REDUCED PEN-UP TIME.
!
!*    IF (ST0) CALL LINE (X0,Y0,X(N1),Y(N1))
!*    IF (.NOT. ST0) CALL LINE (X(N1),Y(N1),X0,Y0)
      if ( n1>=n0 ) st0 = .not.st0
!
!   TEST FOR TERMINATION OF LOOP ON NEIGHBORS
!
      if ( indx==indf ) then
!
!   MARK N0 AS HAVING BEEN PROCESSED
!
         iend(n0) = -indl
         do
!
!   SEARCH THE NEIGHBORS OF N0 FOR AN UNMARKED NODE
!     STARTING WITH IADJ(INDX) = IADJ(INDF)
!
            new = iadj(indx)
            if ( new==0 ) exit
            if ( iend(new)<0 ) then
!
!   TEST FOR TERMINATION
!
               if ( indx==indl ) exit
               indx = indx + 1
            else
               n0 = new
               indf = iabs(iend(n0-1)) + 1
               goto 100
            endif
         enddo
!
!   ALL NEIGHBORS OF N0 ARE MARKED.  SEARCH IEND FOR AN
!     UNMARKED NODE.
!
         do n0 = 2 , nn
            if ( iend(n0)>=0 ) then
               indf = -iend(n0-1) + 1
               goto 100
            endif
         enddo
!
! ALL NODES HAVE BEEN MARKED.  RESTORE IEND.
!
         do n0 = 1 , nn
            iend(n0) = -iend(n0)
         enddo
         if ( numbr==0 ) then
         endif
      else
         indx = indx - 1
         goto 200
      endif
!
! PLOT THE NODAL INDICES
!
!*    THE NUMBR OPTION MAY BE IMPLEMENTED BY INSERTING A
!*    LOOP ON THE NODAL COORDINATES WHICH WRITES INDICES
!*    NEXT TO THE NODES.
!
! TERMINATE PLOTTING -- MOVE TO A NEW FRAME.
!
!*    CALL FRAME
      return
99999 end subroutine trplot
!*==TRPRNT.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine trprnt(n,lunit,x,y,iadj,iend,iflag)
      implicit none
!*--TRPRNT3356
      integer n , lunit , iadj(*) , iend(n) , iflag
      real x(n) , y(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
! THIS ROUTINE PRINTS THE ADJACENCY LISTS AND, OPTIONALLY,
! THE NODAL COORDINATES.  THE NUMBERS OF BOUNDARY NODES,
! TRIANGLES, AND ARCS ARE ALSO PRINTED.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            3 .LE. N .LE. 9999.
!
!                    LUNIT - LOGICAL UNIT FOR OUTPUT.  1
!                            .LE. LUNIT .LE. 99.  OUTPUT IS
!                            PRINTED ON UNIT 6 IF LUNIT IS
!                            OUT OF RANGE.
!
!                      X,Y - VECTORS OF COORDINATES OF THE
!                            NODES IN THE MESH.  NOT USED
!                            UNLESS IFLAG = 0.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF X AND Y ARE TO BE
!                                      PRINTED (TO 6 DECIMAL
!                                      PLACES).
!                            IFLAG = 1 IF X AND Y ARE NOT
!                                      TO BE PRINTED.
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
!               FOR THE LAST PAGE WHICH HAS 3 ADDITIONAL
!               LINES
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
! PRINT X, Y, AND IADJ
!
            write (lun,99003)
99003       format (' ','NODE',5x,'X(NODE)',8x,'Y(NODE)',20x,           &
     &              'NEIGHBORS OF NODE'//)
            nb = 0
            indf = 1
            do node = 1 , nn
               indl = iend(node)
               if ( iadj(indl)==0 ) nb = nb + 1
               inc = (indl-indf)/8 + 2
               nl = nl + inc
               if ( nl>nlmax ) write (lun,99009)
               if ( nl>nlmax ) nl = inc
               write (lun,99004) node , x(node) , y(node) ,             &
     &                           (iadj(i),i=indf,indl)
99004          format (' ',i4,2E15.6,5x,8I5/(' ',39x,8I5))
               if ( indl-indf/=7 ) write (lun,99008)
               indf = indl + 1
            enddo
         else
!
! PRINT IADJ ONLY
!
            write (lun,99005)
99005       format (' ','NODE',32x,'NEIGHBORS OF NODE'//)
            nb = 0
            indf = 1
            do node = 1 , nn
               indl = iend(node)
               if ( iadj(indl)==0 ) nb = nb + 1
               inc = (indl-indf)/14 + 2
               nl = nl + inc
               if ( nl>nlmax ) write (lun,99009)
               if ( nl>nlmax ) nl = inc
               write (lun,99006) node , (iadj(i),i=indf,indl)
99006          format (' ',i4,5x,14I5/(' ',9x,14I5))
               if ( indl-indf/=13 ) write (lun,99008)
               indf = indl + 1
            enddo
         endif
      endif
!
! PRINT NB, NA, AND NT
!
      nt = 2*nn - nb - 2
      na = nt + nn - 1
      write (lun,99007) nb , na , nt
99007 format (/' ','NB = ',i4,' BOUNDARY NODES',10x,'NA = ',i5,' ARCS', &
     &        10x,'NT = ',i5,' TRIANGLES')
      return
99008 format (' ')
99009 format ('1')
      end subroutine trprnt
!*==COORDS.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine coords(x,y,x1,x2,x3,y1,y2,y3,r,ier)
      implicit none
!*--COORDS3506
      integer ier
      real x , y , x1 , x2 , x3 , y1 , y2 , y3 , r(3)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE COMPUTES THE THREE BARYCENTRIC COORDINATES
! OF A POINT IN THE PLANE FOR A GIVEN TRIANGLE.
!
! INPUT PARAMETERS - X,Y - X AND Y COORDINATES OF THE POINT
!                          WHOSE BARYCENTRIC COORDINATES ARE
!                          DESIRED.
!
!      X1,X2,X3,Y1,Y2,Y3 - COORDINATES OF THE VERTICES OF
!                          THE TRIANGLE.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -  R - 3-VECTOR OF BARYCENTRIC COORDI-
!                          NATES UNLESS IER = 1.  NOTE THAT
!                          R(I) .LT. 0. IFF (X,Y) IS TO THE
!                          RIGHT OF THE VECTOR FROM VERTEX
!                          I+1 TO VERTEX I+2 (CYCLICAL
!                          ARITHMETIC).
!
!                    IER - ERROR INDICATOR
!                          IER = 0 IF NO ERRORS WERE
!                                  ENCOUNTERED.
!                          IER = 1 IF THE VERTICES OF THE
!                                  TRIANGLE ARE COLLINEAR.
!
! MODULES REFERENCED BY COORDS - NONE
!
!***********************************************************
!
      real u(3) , v(3) , area , xp , yp
!
! LOCAL PARAMETERS -
!
! U(K),V(K) = X AND Y COMPONENTS OF THE VECTOR REPRESENTING
!               THE SIDE OPPOSITE VERTEX K FOR K = 1,2,3.
! AREA =      TWICE THE AREA OF THE TRIANGLE.
! XP,YP =     X-X1, Y-Y1
!
      u(1) = x3 - x2
      u(2) = x1 - x3
      u(3) = x2 - x1
!
      v(1) = y3 - y2
      v(2) = y1 - y3
      v(3) = y2 - y1
!
! AREA = 3-1 X 3-2
!
      area = u(1)*v(2) - u(2)*v(1)
      if ( area==0. ) then
!
! VERTICES ARE COLLINEAR
!
         ier = 1
         goto 99999
      endif
!
! R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
!   R(3) = (1-2 X 1-(X,Y))/AREA
!
      r(1) = (u(1)*(y-y2)-v(1)*(x-x2))/area
      xp = x - x1
      yp = y - y1
      r(2) = (u(2)*yp-v(2)*xp)/area
      r(3) = (u(3)*yp-v(3)*xp)/area
      ier = 0
      return
99999 end subroutine coords
!*==GIVENS.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine givens(a,b,c,s)
      implicit none
!*--GIVENS3587
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
!*==GRADG.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine gradg(n,x,y,z,iadj,iend,eps,nit,zxzy,ier)
      implicit none
!*--GRADG3685
      integer n , iadj(*) , iend(n) , nit , ier
      real x(n) , y(n) , z(n) , eps , zxzy(2,n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF N NODES IN THE PLANE WITH
! ASSOCIATED DATA VALUES, THIS ROUTINE USES A GLOBAL METHOD
! TO COMPUTE ESTIMATED GRADIENTS AT THE NODES.  THE METHOD
! CONSISTS OF MINIMIZING A QUADRATIC FUNCTIONAL Q(G) OVER
! THE N-VECTOR G OF GRADIENTS WHERE Q APPROXIMATES THE LIN-
! EARIZED CURVATURE OF AN INTERPOLANT F OVER THE TRIANGULA-
! TION.  THE RESTRICTION OF F TO AN ARC OF THE TRIANGULATION
! IS TAKEN TO BE THE HERMITE CUBIC INTERPOLANT OF THE DATA
! VALUES AND TANGENTIAL GRADIENT COMPONENTS AT THE END-
! POINTS OF THE ARC, AND Q IS THE SUM OF THE LINEARIZED
! CURVATURES OF F ALONG THE ARCS -- THE INTEGRALS OVER THE
! ARCS OF D2F(T)**2 WHERE D2F(T) IS THE SECOND DERIVATIVE
! OF F WITH RESPECT TO DISTANCE T ALONG THE ARC.  THIS MIN-
! IMIZATION PROBLEM CORRESPONDS TO AN ORDER 2N SYMMETRIC
! POSITIVE-DEFINITE SPARSE LINEAR SYSTEM WHICH IS SOLVED FOR
! THE X AND Y PARTIAL DERIVATIVES BY THE BLOCK GAUSS-SEIDEL
! METHOD WITH 2 BY 2 BLOCKS.
!   AN ALTERNATIVE METHOD, SUBROUTINE GRADL, COMPUTES A
! LOCAL APPROXIMATION TO THE PARTIALS AT A SINGLE NODE AND
! MAY BE MORE ACCURATE, DEPENDING ON THE DATA VALUES AND
! DISTRIBUTION OF NODES (NEITHER METHOD EMERGED AS SUPERIOR
! IN TESTS FOR ACCURACY).  HOWEVER, IN TESTS RUN ON AN IBM
! 370, GRADG WAS FOUND TO BE ABOUT 3.6 TIMES AS FAST FOR
! NIT = 4.
!
! INPUT PARAMETERS - N - NUMBER OF NODES.  N .GE. 3.
!
!                  X,Y - CARTESIAN COORDINATES OF THE NODES.
!
!                    Z - DATA VALUES AT THE NODES.  Z(I) IS
!                        ASSOCIATED WITH (X(I),Y(I)).
!
!            IADJ,IEND - DATA STRUCTURE DEFINING THE TRIAN-
!                        GULATION.  SEE SUBROUTINE TRMESH.
!
!                  EPS - NONNEGATIVE CONVERGENCE CRITERION.
!                        THE METHOD IS TERMINATED WHEN THE
!                        MAXIMUM CHANGE IN A GRADIENT COMPO-
!                        NENT BETWEEN ITERATIONS IS AT MOST
!                        EPS.  EPS = 1.E-2 IS SUFFICIENT FOR
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
!                        AND NIT = 3 OR 4.
!
!                 ZXZY - 2 BY N ARRAY WHOSE COLUMNS CONTAIN
!                        INITIAL ESTIMATES OF THE PARTIAL
!                        DERIVATIVES (ZERO VECTORS ARE
!                        SUFFICIENT).
!
! OUTPUT PARAMETERS - NIT - NUMBER OF GAUSS-SEIDEL ITERA-
!                           TIONS EMPLOYED.
!
!                    ZXZY - ESTIMATED X AND Y PARTIAL DERIV-
!                           ATIVES AT THE NODES WITH X PAR-
!                           TIALS IN THE FIRST ROW.  ZXZY IS
!                           NOT CHANGED IF IER = 2.
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
! MODULES REFERENCED BY GRADG - NONE
!
! INTRINSIC FUNCTIONS CALLED BY GRADG - SQRT, AMAX1, ABS
!
!***********************************************************
!
      integer nn , maxit , iter , k , indf , indl , indx , nb
      real tol , dgmax , xk , yk , zk , zxk , zyk , a11 , a12 , a22 ,   &
     &     r1 , r2 , delx , dely , delxs , delys , dsq , dcub , t ,     &
     &     dzx , dzy
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
! XK,YK,ZK =    X(K), Y(K), Z(K)
! ZXK,ZYK =     INITIAL VALUES OF ZXZY(1,K) AND ZXZY(2,K)
! A11,A12,A22 = MATRIX COMPONENTS OF THE 2 BY 2 BLOCK A*DG
!                 = R WHERE A IS SYMMETRIC, DG = (DZX,DZY)
!                 IS THE CHANGE IN THE GRADIENT AT K, AND R
!                 IS THE RESIDUAL
! R1,R2 =       COMPONENTS OF THE RESIDUAL -- DERIVATIVES OF
!                 Q WITH RESPECT TO THE COMPONENTS OF THE
!                 GRADIENT AT NODE K
! DELX,DELY =   COMPONENTS OF THE ARC NB-K
! DELXS,DELYS = DELX**2, DELY**2
! DSQ =         SQUARE OF THE DISTANCE D BETWEEN K AND NB
! DCUB =        D**3
! T =           FACTOR OF R1 AND R2
! DZX,DZY =     SOLUTION OF THE 2 BY 2 SYSTEM -- CHANGE IN
!                 DERIVATIVES AT K FROM THE PREVIOUS ITERATE
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
               zxk = zxzy(1,k)
               zyk = zxzy(2,k)
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
!   COMPUTE THE COMPONENTS OF ARC NB-K
!
                     delx = x(nb) - xk
                     dely = y(nb) - yk
                     delxs = delx*delx
                     delys = dely*dely
                     dsq = delxs + delys
                     dcub = dsq*sqrt(dsq)
!
!   UPDATE THE SYSTEM COMPONENTS FOR NODE NB
!
                     a11 = a11 + delxs/dcub
                     a12 = a12 + delx*dely/dcub
                     a22 = a22 + delys/dcub
                     t = (1.5*(z(nb)-zk)                                &
     &                   -((zxzy(1,nb)/2.+zxk)*delx+(zxzy(2,nb)/2.+zyk) &
     &                   *dely))/dcub
                     r1 = r1 + t*delx
                     r2 = r2 + t*dely
                  endif
               enddo
!
!   SOLVE THE 2 BY 2 SYSTEM AND UPDATE DGMAX
!
               dzy = (a11*r2-a12*r1)/(a11*a22-a12*a12)
               dzx = (r1-a12*dzy)/a11
               dgmax = amax1(dgmax,abs(dzx),abs(dzy))
!
!   UPDATE THE PARTIALS AT NODE K
!
               zxzy(1,k) = zxk + dzx
               zxzy(2,k) = zyk + dzy
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
!*==GRADL.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine gradl(n,k,x,y,z,iadj,iend,dx,dy,ier)
      implicit none
!*--GRADL3911
      integer n , k , iadj(1) , iend(n) , ier
      real x(n) , y(n) , z(n) , dx , dy
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF N POINTS IN THE PLANE
! WITH ASSOCIATED DATA VALUES Z, THIS SUBROUTINE ESTIMATES
! X AND Y PARTIAL DERIVATIVES AT NODE K.  THE DERIVATIVES
! ARE TAKEN TO BE THE PARTIALS AT K OF A QUADRATIC FUNCTION
! WHICH INTERPOLATES Z(K) AND FITS THE DATA VALUES AT A SET
! OF NEARBY NODES IN A WEIGHTED LEAST SQUARES SENSE. A MAR-
! QUARDT STABILIZATION FACTOR IS USED IF NECESSARY TO ENSURE
! A WELL-CONDITIONED SYSTEM AND A LINEAR FITTING FUNCTION IS
! USED IF N .LT. 6.  THUS, A UNIQUE SOLUTION EXISTS UNLESS
! THE NODES ARE COLLINEAR.
!   AN ALTERNATIVE ROUTINE, GRADG, EMPLOYS A GLOBAL METHOD
! TO COMPUTE THE PARTIAL DERIVATIVES AT ALL OF THE NODES AT
! ONCE.  THAT METHOD IS MORE EFFICIENT (WHEN ALL PARTIALS
! ARE NEEDED) AND MAY BE MORE ACCURATE, DEPENDING ON THE
! DATA.
!
! INPUT PARAMETERS - N - NUMBER OF NODES IN THE TRIANGULA-
!                        TION.  N .GE. 3.
!
!                    K - NODE AT WHICH DERIVATIVES ARE
!                        SOUGHT.  1 .LE. K .LE. N.
!
!                  X,Y - N-VECTORS CONTAINING THE CARTESIAN
!                        COORDINATES OF THE NODES.
!
!                    Z - N-VECTOR CONTAINING THE DATA VALUES
!                        ASSOCIATED WITH THE NODES.
!
!                 IADJ - SET OF ADJACENCY LISTS.
!
!                 IEND - POINTERS TO THE ENDS OF ADJACENCY
!                        LISTS FOR EACH NODE.
!
! IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS - DX,DY - ESTIMATED PARTIAL DERIVATIVES
!                             AT NODE K UNLESS IER .LT. 0.
!
!                       IER - ERROR INDICATOR
!                             IER .GT. 0 IF NO ERRORS WERE
!                                      ENCOUNTERED.  IER
!                                      CONTAINS THE NUMBER
!                                      OF NODES (INCLUDING
!                                      K) USED IN THE FIT.
!                                      IER = 3, 4, OR 5 IM-
!                                      PLIES A LINEAR FIT.
!                             IER = -1 IF N OR K IS OUT OF
!                                      RANGE.
!                             IER = -2 IF ALL NODES ARE
!                                      COLLINEAR.
!
! MODULES REFERENCED BY GRADL - GETNP, SETUP, GIVENS,
!                               ROTATE
!
! INTRINSIC FUNCTIONS CALLED BY GRADL - MIN0, FLOAT, SQRT,
!                                       AMIN1, ABS
!
!***********************************************************
!
      integer nn , kk , lmn , lmx , lmin , lmax , lm1 , lnp , npts(30) ,&
     &        ierr , np , i , j , im1 , jp1 , ip1 , l
      real sum , ds , r , rs , rtol , avsq , av , xk , yk , zk , a(6,6) &
     &     , c , s , dmin , dtol , sf
      data lmn/10/
      data lmx/30/ , rtol/1.E-5/ , dtol/.01/ , sf/1./
!
! LOCAL PARAMETERS -
!
! NN,KK =     LOCAL COPIES OF N AND K
! LMN,LMX =   MINIMUM AND MAXIMUM VALUES OF LNP FOR N
!               SUFFICIENTLY LARGE.  IN MOST CASES LMN-1
!               NODES ARE USED IN THE FIT.  4 .LE. LMN .LE.
!               LMX.
! LMIN,LMAX = MIN(LMN,N), MIN(LMX,N)
! LM1 =       LMIN-1 OR LNP-1
! LNP =       LENGTH OF NPTS
! NPTS =      ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!               NODES ORDERED BY DISTANCE FROM K.  NPTS(1)=K
!               AND THE FIRST LNP-1 ELEMENTS OF NPTS ARE
!               USED IN THE LEAST SQUARES FIT.  UNLESS LNP
!               EXCEEDS LMAX, NPTS(LNP) DETERMINES R.
! IERR =      ERROR FLAG FOR CALLS TO GETNP (NOT CHECKED)
! NP =        ELEMENT OF NPTS TO BE ADDED TO THE SYSTEM
! I,J =       DO-LOOP INDICES
! IM1,JP1 =   I-1, J+1
! IP1 =       I+1
! L =         NUMBER OF COLUMNS OF A**T TO WHICH A ROTATION
!               IS APPLIED
! SUM =       SUM OF SQUARED EUCLIDEAN DISTANCES BETWEEN
!               NODE K AND THE NODES USED IN THE LEAST
!               SQUARES FIT
! DS =        SQUARED DISTANCE BETWEEN NODE K AND AN ELE-
!               MENT OF NPTS
! R =         DISTANCE BETWEEN NODE K AND NPTS(LNP) OR SOME
!               POINT FURTHER FROM K THAN NPTS(LMAX) IF
!               NPTS(LMAX) IS USED IN THE FIT.  R IS A
!               RADIUS OF INFLUENCE WHICH ENTERS INTO THE
!               WEIGHTS (SEE SUBROUTINE SETUP).
! RS =        R*R
! RTOL =      TOLERANCE FOR DETERMINING R.  IF THE RELATIVE
!               CHANGE IN DS BETWEEN TWO ELEMENTS OF NPTS IS
!               NOT GREATER THAN RTOL THEY ARE TREATED AS
!               BEING THE SAME DISTANCE FROM NODE K
! AVSQ =      AV*AV
! AV =        ROOT-MEAN-SQUARE DISTANCE BETWEEN K AND THE
!               NODES (OTHER THAN K) IN THE LEAST SQUARES
!               FIT.  THE FIRST 3 COLUMNS OF THE SYSTEM ARE
!               SCALED BY 1/AVSQ, THE NEXT 2 BY 1/AV.
! XK,YK,ZK =  COORDINATES AND DATA VALUE ASSOCIATED WITH K
! A =         TRANSPOSE OF THE AUGMENTED REGRESSION MATRIX
! C,S =       COMPONENTS OF THE PLANE ROTATION DETERMINED
!               BY SUBROUTINE GIVENS
! DMIN =      MINIMUM OF THE MAGNITUDES OF THE DIAGONAL
!               ELEMENTS OF THE REGRESSION MATRIX AFTER
!               ZEROS ARE INTRODUCED BELOW THE DIAGONAL
! DTOL =      TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!               SYSTEM.  THE SYSTEM IS ACCEPTED WHEN DMIN
!               .GE. DTOL
! SF =        MARQUARDT STABILIZATION FACTOR USED TO DAMP
!               OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!               PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!               IS ILL-CONDITIONED.  AS SF INCREASES, THE
!               FITTING FUNCTION APPROACHES A LINEAR
!
      nn = n
      kk = k
!
! CHECK FOR ERRORS AND INITIALIZE LMIN, LMAX
!
      if ( nn<3 .or. kk<1 .or. kk>nn ) then
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
!   SET NPTS TO THE CLOSEST LMIN-1 NODES TO K.
!
         sum = 0.
         npts(1) = kk
         lm1 = lmin - 1
         do lnp = 2 , lm1
            call getnp(x,y,iadj,iend,lnp,npts,ds,ierr)
            sum = sum + ds
         enddo
!
! ADD ADDITIONAL NODES TO NPTS UNTIL THE RELATIVE INCREASE
!   IN DS IS AT LEAST RTOL.
!
         do lnp = lmin , lmax
            call getnp(x,y,iadj,iend,lnp,npts,rs,ierr)
            if ( (rs-ds)/ds>rtol ) then
               if ( lnp>6 ) goto 50
            endif
            sum = sum + rs
         enddo
!
! USE ALL LMAX NODES IN THE LEAST SQUARES FIT.  RS IS
!   ARBITRARILY INCREASED BY 10 PER CENT.
!
         rs = 1.1*rs
         lnp = lmax + 1
!
! THERE ARE LNP-2 EQUATIONS CORRESPONDING TO NODES NPTS(2),
!   ...,NPTS(LNP-1).
!
 50      avsq = sum/float(lnp-2)
         av = sqrt(avsq)
         r = sqrt(rs)
         xk = x(kk)
         yk = y(kk)
         zk = z(kk)
         if ( lnp<7 ) then
!
! 4 .LE. LNP .LE. 6 (2, 3, OR 4 EQUATIONS) -- FIT A PLANE TO
!   THE DATA USING THE LAST 3 COLUMNS OF A.
!
            np = npts(2)
            call setup(xk,yk,zk,x(np),y(np),z(np),av,avsq,r,a(1,4))
            np = npts(3)
            call setup(xk,yk,zk,x(np),y(np),z(np),av,avsq,r,a(1,5))
            call givens(a(4,4),a(4,5),c,s)
            call rotate(2,c,s,a(5,4),a(5,5))
            if ( lnp/=4 ) then
!
               lm1 = lnp - 1
               do i = 4 , lm1
                  np = npts(i)
                  call setup(xk,yk,zk,x(np),y(np),z(np),av,avsq,r,a(1,6)&
     &                       )
                  call givens(a(4,4),a(4,6),c,s)
                  call rotate(2,c,s,a(5,4),a(5,6))
                  call givens(a(5,5),a(5,6),c,s)
                  call rotate(1,c,s,a(6,5),a(6,6))
               enddo
            endif
         else
!
! SET UP THE FIRST 5 EQUATIONS OF THE AUGMENTED REGRESSION
!   MATRIX (TRANSPOSED) AS THE COLUMNS OF A, AND ZERO OUT
!   THE LOWER TRIANGLE (UPPER TRIANGLE OF A) WITH GIVENS
!   ROTATIONS
!
            do i = 1 , 5
               np = npts(i+1)
               call setup(xk,yk,zk,x(np),y(np),z(np),av,avsq,r,a(1,i))
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
                  dmin = amin1(abs(a(1,1)),abs(a(2,2)),abs(a(3,3)),     &
     &                   abs(a(4,4)),abs(a(5,5)))
                  if ( dmin>=dtol ) goto 100
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
                     exit
                  else
!
! ADD ANOTHER NODE TO THE SYSTEM AND INCREASE R --
!   I .EQ. LNP
!
                     lnp = lnp + 1
                     if ( lnp<=lmax )                                   &
     &                    call getnp(x,y,iadj,iend,lnp,npts,rs,ierr)
                     r = sqrt(1.1*rs)
                  endif
               else
                  np = npts(i)
                  call setup(xk,yk,zk,x(np),y(np),z(np),av,avsq,r,a(1,6)&
     &                       )
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
! TEST THE LINEAR FIT FOR ILL-CONDITIONING
!
         dmin = amin1(abs(a(4,4)),abs(a(5,5)))
         if ( dmin<dtol ) then
!
! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES
!
            ier = -2
            goto 99999
         endif
      endif
!
! SOLVE THE 2 BY 2 TRIANGULAR SYSTEM FOR THE DERIVATIVES
!
 100  dy = a(6,5)/a(5,5)
      dx = (a(6,4)-a(5,4)*dy)/a(4,4)/av
      dy = dy/av
      ier = lnp - 1
      return
99999 end subroutine gradl
!*==INTRC0.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine intrc0(n,px,py,x,y,z,iadj,iend,ist,pz,ier)
      implicit none
!*--INTRC04224
      integer n , iadj(*) , iend(n) , ist , ier
      real px , py , x(n) , y(n) , z(n) , pz
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
! THIS ROUTINE COMPUTES THE VALUE AT (PX,PY) OF A PIECEWISE
! LINEAR SURFACE WHICH INTERPOLATES DATA VALUES AT THE
! VERTICES OF THE TRIANGLES.  THE SURFACE IS EXTENDED IN A
! CONTINUOUS FASHION BEYOND THE BOUNDARY OF THE TRIANGULAR
! MESH, ALLOWING EXTRAPOLATION.  INTRC0 IS PART OF AN
! INTERPOLATION PACKAGE WHICH PROVIDES ROUTINES TO GENERATE,
! UPDATE, AND PLOT THE MESH.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            N .GE. 3.
!
!                    PX,PY - POINT AT WHICH THE INTERPOLATED
!                            VALUE IS DESIRED.
!
!                      X,Y - VECTORS OF COORDINATES OF THE
!                            NODES IN THE MESH.
!
!                        Z - VECTOR OF DATA VALUES AT THE
!                            NODES.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING (PX,PY).  1 .LE. IST
!                            .LE. N.  THE OUTPUT VALUE OF
!                            IST FROM A PREVIOUS CALL MAY
!                            BE A GOOD CHOICE.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS -  IST - INDEX OF ONE OF THE VERTICES OF
!                            THE TRIANGLE CONTAINING (PX,PY)
!                            UNLESS IER .LT. 0.
!
!                       PZ - VALUE OF THE INTERPOLATORY
!                            SURFACE AT (PX,PY) OR ZERO
!                            IF IER .LT. 0.
!
!                      IER - ERROR INDICATOR
!                            IER = 0 IF NO ERRORS WERE
!                                    ENCOUNTERED.
!                            IER = 1 IF NO ERRORS WERE EN-
!                                    COUNTERED AND EXTRAPO-
!                                    LATION WAS PERFORMED.
!                            IER = -1 IF N OR IST IS OUT OF
!                                     RANGE.
!                            IER = -2 IF THE NODES ARE COL-
!                                     LINEAR.
!
! MODULES REFERENCED BY INTRC0 - TRFIND, COORDS
!
!***********************************************************
!
      integer i1 , i2 , i3 , n1 , n2 , indx
      real xp , yp , r(3) , x1 , y1 , x2 , y2 , dp
!
! LOCAL PARAMETERS -
!
! I1,I2,I3 = VERTEX INDICES RETURNED BY TRFIND
! N1,N2 =    ENDPOINTS OF THE CLOSEST BOUNDARY EDGE TO P
!              WHEN P IS OUTSIDE OF THE MESH BOUNDARY
! INDX =     IADJ INDEX OF N1 AS A NEIGHBOR OF N2
! XP,YP =    LOCAL COPIES OF THE COORDINATES OF P=(PX,PY)
! R =        BARYCENTRIC COORDINATES
! X1,Y1 =    X,Y COORDINATES OF N1
! X2,Y2 =    X,Y COORDINATES OF N2
! DP =       INNER PRODUCT OF N1-N2 AND P-N2
!
      if ( n<3 .or. ist<1 .or. ist>n ) then
!
! N .LT. 3 OR IST IS OUT OF RANGE
!
         pz = 0.
         ier = -1
         return
      else
         xp = px
         yp = py
!
! FIND A TRIANGLE CONTAINING P IF P IS WITHIN THE MESH
!   BOUNDARY
!
         call trfind(ist,xp,yp,x,y,iadj,iend,i1,i2,i3)
         if ( i1/=0 ) then
            ist = i1
            if ( i3==0 ) then
!
! P IS OUTSIDE OF THE MESH BOUNDARY.  EXTRAPOLATE TO P BY
!   EXTENDING THE INTERPOLATORY SURFACE AS A CONSTANT
!   BEYOND THE BOUNDARY.  THUS PZ IS THE SURFACE FUNCTION
!   VALUE AT Q WHERE Q IS THE CLOSEST BOUNDARY POINT TO P.
!
! DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM THE
!   RIGHTMOST VISIBLE NODE I1.
!
               n2 = i1
               do
!
! SET N1 TO THE LAST NONZERO NEIGHBOR OF N2 AND COMPUTE DP
!
                  indx = iend(n2) - 1
                  n1 = iadj(indx)
                  x1 = x(n1)
                  y1 = y(n1)
                  x2 = x(n2)
                  y2 = y(n2)
                  dp = (x1-x2)*(xp-x2) + (y1-y2)*(yp-y2)
                  if ( dp<=0. ) then
!
! N2 IS THE CLOSEST BOUNDARY POINT TO P
!
                     pz = z(n2)
                     ier = 1
                     return
                  elseif ( (xp-x1)*(x2-x1)+(yp-y1)*(y2-y1)>0. ) then
!
! THE CLOSEST BOUNDARY POINT TO P LIES ON N2-N1.  COMPUTE
!   ITS COORDINATES WITH RESPECT TO N2-N1.
!
                     r(1) = dp/((x2-x1)**2+(y2-y1)**2)
                     r(2) = 1. - r(1)
                     pz = r(1)*z(n1) + r(2)*z(n2)
                     ier = 1
                     return
                  else
                     n2 = n1
                  endif
               enddo
            else
!
! COMPUTE BARYCENTRIC COORDINATES
!
               call coords(xp,yp,x(i1),x(i2),x(i3),y(i1),y(i2),y(i3),r, &
     &                     ier)
               if ( ier==0 ) then
                  pz = r(1)*z(i1) + r(2)*z(i2) + r(3)*z(i3)
                  return
               endif
            endif
         endif
!
! NODES ARE COLLINEAR
!
         pz = 0.
         ier = -2
      endif
      end subroutine intrc0
!*==INTRC1.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine intrc1(n,px,py,x,y,z,iadj,iend,iflag,zxzy,ist,pz,ier)
      implicit none
!*--INTRC14394
      integer n , iadj(*) , iend(n) , iflag , ist , ier
      real px , py , x(n) , y(n) , z(n) , zxzy(2,n) , pz
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,
! THIS ROUTINE DETERMINES A PIECEWISE CUBIC FUNCTION F(X,Y)
! WHICH INTERPOLATES A SET OF DATA VALUES AND PARTIAL
! DERIVATIVES AT THE VERTICES.  F HAS CONTINUOUS FIRST
! DERIVATIVES OVER THE MESH AND EXTENDS BEYOND THE MESH
! BOUNDARY ALLOWING EXTRAPOLATION.  INTERPOLATION IS EXACT
! FOR QUADRATIC DATA.  THE VALUE OF F AT (PX,PY) IS
! RETURNED.  INTRC1 IS PART OF AN INTERPOLATION PACKAGE
! WHICH PROVIDES ROUTINES TO GENERATE, UPDATE AND PLOT THE
! MESH.
!
! INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.
!                            N .GE. 3.
!
!                    PX,PY - COORDINATES OF A POINT AT WHICH
!                            F IS TO BE EVALUATED.
!
!                      X,Y - VECTORS OF COORDINATES OF THE
!                            NODES IN THE MESH.
!
!                        Z - VECTOR OF DATA VALUES AT THE
!                            NODES.
!
!                     IADJ - SET OF ADJACENCY LISTS OF NODES
!                            IN THE MESH.
!
!                     IEND - POINTERS TO THE ENDS OF
!                            ADJACENCY LISTS IN IADJ FOR
!                            EACH NODE IN THE MESH.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF INTRC1 IS TO
!                                      PROVIDE DERIVATIVE
!                                      ESTIMATES (FROM
!                                      GRADL).
!                            IFLAG = 1 IF DERIVATIVES ARE
!                                      USER PROVIDED.
!
!                     ZXZY - 2 BY N ARRAY WHOSE COLUMNS
!                            CONTAIN ESTIMATED PARTIAL DER-
!                            IVATIVES AT THE NODES (X PAR-
!                            TIALS IN THE FIRST ROW) IF
!                            IFLAG = 1, NOT USED IF IFLAG
!                            = 0.
!
!                      IST - INDEX OF THE STARTING NODE IN
!                            THE SEARCH FOR A TRIANGLE CON-
!                            TAINING (PX,PY).  1 .LE. IST
!                            .LE. N.  THE OUTPUT VALUE OF
!                            IST FROM A PREVIOUS CALL MAY
!                            BE A GOOD CHOICE.
!
! IADJ AND IEND MAY BE CREATED BY TRMESH AND DERIVATIVE
!   ESTIMATES MAY BE COMPUTED BY GRADL OR GRADG.
!
! INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS
!   ROUTINE.
!
! OUTPUT PARAMETERS - IST - INDEX OF ONE OF THE VERTICES OF
!                           THE TRIANGLE CONTAINING (PX,PY)
!                           UNLESS IER .LT. 0.
!
!                      PZ - VALUE OF F AT (PX,PY), OR 0 IF
!                           IER .LT. 0.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE
!                                   ENCOUNTERED.
!                           IER = 1 IF NO ERRORS WERE EN-
!                                   COUNTERED AND EXTRAPOLA-
!                                   TION WAS PERFORMED.
!                           IER = -1 IF N, IFLAG, OR IST IS
!                                    OUT OF RANGE.
!                           IER = -2 IF THE NODES ARE COL-
!                                    LINEAR.
!
! MODULES REFERENCED BY INTRC1 - TRFIND, TVAL,
!             (AND OPTIONALLY)   GRADL, GETNP, SETUP,
!                                GIVENS, ROTATE
!
!***********************************************************
!
      integer nn , i1 , i2 , i3 , ierr , n1 , n2 , indx
      real xp , yp , zx1 , zy1 , zx2 , zy2 , zx3 , zy3 , x1 , y1 , x2 , &
     &     y2 , x3 , y3 , z1 , z2 , z3 , dum , dp , u , v , xq , yq ,   &
     &     r1 , r2 , a1 , a2 , b1 , b2 , c1 , c2 , f1 , f2
!
! LOCAL PARAMETERS -
!
! NN =                      LOCAL COPY OF N
! I1,I2,I3 =                VERTICES DETERMINED BY TRFIND
! IERR =                    ERROR FLAG FOR CALLS TO GRADL
!                             AND TVAL
! N1,N2 =                   ENDPOINTS OF THE CLOSEST BOUND-
!                             ARY EDGE TO P WHEN P IS OUT-
!                             SIDE OF THE MESH BOUNDARY
! INDX =                    IADJ INDEX OF N1 AS A NEIGHBOR
!                             OF N2
! XP,YP =                   LOCAL COPIES OF THE COORDINATES
!                             OF P=(PX,PY)
! ZX1,ZY1,ZX2,ZY2,ZX3,ZY3 = X AND Y DERIVATIVES AT THE
!                             VERTICES OF A TRIANGLE T WHICH
!                             CONTAINS P OR AT N1 AND N2
! X1,Y1,X2,Y2,X3,Y3 =       X,Y COORDINATES OF THE VERTICES
!                             OF T OR OF N1 AND N2
! Z1,Z2,Z3 =                DATA VALUES AT THE VERTICES OF T
! DUM =                     DUMMY VARIABLE FOR CALL TO TVAL
! DP =                      INNER PRODUCT OF N1-N2 AND P-N2
! U,V =                     X,Y COORDINATES OF THE VECTOR
!                             N2-N1
! XQ,YQ =                   X,Y COORDINATES OF THE CLOSEST
!                             BOUNDARY POINT TO P WHEN P IS
!                             OUTSIDE OF THE MESH BOUNDARY
! R1,R2 =                   BARYCENTRIC COORDINATES OF Q
!                             WITH RESPECT TO THE LINE SEG-
!                             MENT N2-N1 CONTAINING Q
! A1,A2,B1,B2,C1,C2 =       CARDINAL FUNCTIONS FOR EVALUAT-
!                             ING THE INTERPOLATORY SURFACE
!                             AT Q
! F1,F2 =                   CUBIC FACTORS USED TO COMPUTE
!                             THE CARDINAL FUNCTIONS
!
      nn = n
      pz = 0.
      if ( nn<3 .or. iflag<0 .or. iflag>1 .or. ist<1 .or. ist>nn ) then
!
! N, IFLAG, OR IST OUT OF RANGE
!
         ier = -1
         return
      else
         xp = px
         yp = py
!
! FIND A TRIANGLE CONTAINING P IF P IS WITHIN THE MESH
!   BOUNDARY
!
         call trfind(ist,xp,yp,x,y,iadj,iend,i1,i2,i3)
         if ( i1==0 ) then
!
! NODES ARE COLLINEAR
!
            ier = -2
         else
            ist = i1
            if ( i3==0 ) then
!
! P IS OUTSIDE OF THE MESH BOUNDARY.  EXTRAPOLATE TO P BY
!   PASSING A LINEAR FUNCTION OF ONE VARIABLE THROUGH THE
!   VALUE AND DIRECTIONAL DERIVATIVE (IN THE DIRECTION
!   P-Q) OF THE INTERPOLATORY SURFACE (TVAL) AT Q WHERE
!   Q IS THE CLOSEST BOUNDARY POINT TO P.
!
! DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM
!   THE RIGHTMOST VISIBLE NODE I1.
!
               n2 = i1
               do
!
! SET N1 TO THE LAST NONZERO NEIGHBOR OF N2 AND COMPUTE DP
!
                  indx = iend(n2) - 1
                  n1 = iadj(indx)
                  x1 = x(n1)
                  y1 = y(n1)
                  x2 = x(n2)
                  y2 = y(n2)
                  dp = (x1-x2)*(xp-x2) + (y1-y2)*(yp-y2)
                  if ( dp<=0. ) then
!
! N2 IS THE CLOSEST BOUNDARY POINT TO P.  COMPUTE PARTIAL
!   DERIVATIVES AT N2.
!
                     if ( iflag/=1 ) then
                        call gradl(nn,n2,x,y,z,iadj,iend,zx2,zy2,ierr)
                     else
                        zx2 = zxzy(1,n2)
                        zy2 = zxzy(2,n2)
                     endif
!
! COMPUTE EXTRAPOLATED VALUE AT P
!
                     pz = z(n2) + zx2*(xp-x2) + zy2*(yp-y2)
                     ier = 1
                     return
                  elseif ( (xp-x1)*(x2-x1)+(yp-y1)*(y2-y1)>0. ) then
!
! THE CLOSEST BOUNDARY POINT Q LIES ON N2-N1.  COMPUTE
!   PARTIALS AT N1 AND N2.
!
                     if ( iflag/=1 ) then
                        call gradl(nn,n1,x,y,z,iadj,iend,zx1,zy1,ierr)
                        call gradl(nn,n2,x,y,z,iadj,iend,zx2,zy2,ierr)
                     else
                        zx1 = zxzy(1,n1)
                        zy1 = zxzy(2,n1)
                        zx2 = zxzy(1,n2)
                        zy2 = zxzy(2,n2)
                     endif
!
! COMPUTE Q, ITS BARYCENTRIC COORDINATES, AND THE CARDINAL
!   FUNCTIONS FOR EXTRAPOLATION
!
                     u = x2 - x1
                     v = y2 - y1
                     r1 = dp/(u**2+v**2)
                     r2 = 1. - r1
                     xq = r1*x1 + r2*x2
                     yq = r1*y1 + r2*y2
                     f1 = r1*r1*r2
                     f2 = r1*r2*r2
                     a1 = r1 + (f1-f2)
                     a2 = r2 - (f1-f2)
                     b1 = u*f1
                     b2 = -u*f2
                     c1 = v*f1
                     c2 = -v*f2
!
! COMPUTE THE VALUE OF THE INTERPOLATORY SURFACE (TVAL)
!   AT Q
!
                     pz = a1*z(n1) + a2*z(n2) + b1*zx1 + b2*zx2 +       &
     &                    c1*zy1 + c2*zy2
!
! COMPUTE THE EXTRAPOLATED VALUE AT P
!
                     pz = pz + (r1*zx1+r2*zx2)*(xp-xq) + (r1*zy1+r2*zy2)&
     &                    *(yp-yq)
                     ier = 1
                     return
                  else
                     n2 = n1
                  endif
               enddo
            else
               if ( iflag/=1 ) then
!
! COMPUTE DERIVATIVE ESTIMATES AT THE VERTICES
!
                  call gradl(nn,i1,x,y,z,iadj,iend,zx1,zy1,ierr)
                  call gradl(nn,i2,x,y,z,iadj,iend,zx2,zy2,ierr)
                  call gradl(nn,i3,x,y,z,iadj,iend,zx3,zy3,ierr)
               else
!
! DERIVATIVES ARE USER PROVIDED
!
                  zx1 = zxzy(1,i1)
                  zx2 = zxzy(1,i2)
                  zx3 = zxzy(1,i3)
                  zy1 = zxzy(2,i1)
                  zy2 = zxzy(2,i2)
                  zy3 = zxzy(2,i3)
               endif
!
! SET LOCAL PARAMETERS FOR CALL TO TVAL
!
               x1 = x(i1)
               y1 = y(i1)
               x2 = x(i2)
               y2 = y(i2)
               x3 = x(i3)
               y3 = y(i3)
               z1 = z(i1)
               z2 = z(i2)
               z3 = z(i3)
               call tval(xp,yp,x1,x2,x3,y1,y2,y3,z1,z2,z3,zx1,zx2,zx3,  &
     &                   zy1,zy2,zy3,0,pz,dum,dum,ierr)
               if ( ierr/=0 ) then
                  ier = -2
               else
                  ier = 0
                  return
               endif
            endif
         endif
      endif
      end subroutine intrc1
!*==ROTATE.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine rotate(n,c,s,x,y)
      implicit none
!*--ROTATE4684
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
! LINPACK.
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
!*==SETUP.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine setup(xk,yk,zk,xi,yi,zi,s1,s2,r,row)
      implicit none
!*--SETUP4740
      real xk , yk , zk , xi , yi , zi , s1 , s2 , r , row(6)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS ROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED RE-
! GRESSION MATRIX FOR A WEIGHTED LEAST-SQUARES FIT OF A
! QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES Z WHERE
! Q(XK,YK) = ZK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS) ARE
! SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LINEAR
! TERMS) ARE SCALED BY 1/S1.  THE WEIGHT IS (R-D)/(R*D) IF
! R .GT. D AND 0 IF R .LE. D, WHERE D IS THE DISTANCE
! BETWEEN NODES I AND K.
!
! INPUT PARAMETERS - XK,YK,ZK - COORDINATES AND DATA VALUE
!                               AT NODE K -- INTERPOLATED
!                               BY Q.
!
!                    XI,YI,ZI - COORDINATES AND DATA VALUE
!                               AT NODE I.
!
!                       S1,S2 - INVERSE SCALE FACTORS.
!
!                           R - RADIUS OF INFLUENCE ABOUT
!                               NODE K DEFINING THE WEIGHT.
!
!                         ROW - VECTOR OF LENGTH 6.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETER - ROW - VECTOR CONTAINING A ROW OF THE
!                          AUGMENTED REGRESSION MATRIX.
!
! MODULES REFERENCED BY SETUP - NONE
!
! INTRINSIC FUNCTION CALLED BY SETUP - SQRT
!
!***********************************************************
!
      integer i
      real dx , dy , dxsq , dysq , d , w , w1 , w2
!
! LOCAL PARAMETERS -
!
! I =    DO-LOOP INDEX
! DX =   XI - XK
! DY =   YI - YK
! DXSQ = DX*DX
! DYSQ = DY*DY
! D =    DISTANCE BETWEEN NODES K AND I
! W =    WEIGHT ASSOCIATED WITH THE ROW
! W1 =   W/S1
! W2 =   W/S2
!
      dx = xi - xk
      dy = yi - yk
      dxsq = dx*dx
      dysq = dy*dy
      d = sqrt(dxsq+dysq)
      if ( d<=0. .or. d>=r ) then
!
! NODES K AND I COINCIDE OR NODE I IS OUTSIDE OF THE RADIUS
!   OF INFLUENCE.  SET ROW TO THE ZERO VECTOR.
!
         do i = 1 , 6
            row(i) = 0.
         enddo
         goto 99999
      endif
      w = (r-d)/r/d
      w1 = w/s1
      w2 = w/s2
      row(1) = dxsq*w2
      row(2) = dx*dy*w2
      row(3) = dysq*w2
      row(4) = dx*w1
      row(5) = dy*w1
      row(6) = (zi-zk)*w
      return
99999 end subroutine setup
!*==TRVOL.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      function trvol(x1,x2,x3,y1,y2,y3,z1,z2,z3)
      implicit none
!*--TRVOL4827
!*** Start of declarations inserted by SPAG
      real area , t1 , t2 , t3 , trvol
!*** End of declarations inserted by SPAG
      real x1 , x2 , x3 , y1 , y2 , y3 , z1 , z2 , z3
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   THIS FUNCTION COMPUTES THE INTEGRAL OVER A TRIANGLE OF
! THE PLANAR SURFACE WHICH INTERPOLATES DATA VALUES AT THE
! VERTICES.
!
! INPUT PARAMETERS - X1,X2,X3 - X COORDINATES OF THE
!                               VERTICES OF THE TRIANGLE IN
!                               COUNTERCLOCKWISE ORDER.
!
!                    Y1,Y2,Y3 - Y COORDINATES OF THE
!                               VERTICES OF THE TRIANGLE IN
!                               THE SAME ORDER AS THE X
!                               COORDINATES.
!
!                    Z1,Z2,Z3 - DATA VALUES AT THE VERTICES
!                               IN THE SAME ORDER AS THE
!                               COORDINATES.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER -    TRVOL - VOLUME UNDER THE INTERPOLA-
!                               TORY SURFACE ABOVE THE
!                               TRIANGLE OR ZERO IF THE
!                               COORDINATES ARE INCORRECTLY
!                               ORDERED OR COLLINEAR.
!
! MODULES REFERENCED BY TRVOL - NONE
!
!***********************************************************
!
      t1 = x2*y3 - x3*y2
      t2 = x3*y1 - x1*y3
      t3 = x1*y2 - x2*y1
      area = t1 + t2 + t3
      if ( area<0. ) area = 0.
!
! AREA IS TWICE THE AREA OF THE TRIANGLE.  TRVOL IS THE MEAN
!   OF THE DATA VALUES TIMES THE AREA OF THE TRIANGLE.
!
      trvol = (z1+z2+z3)*area/6.
      end function trvol
!*==TVAL.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine tval(x,y,x1,x2,x3,y1,y2,y3,z1,z2,z3,zx1,zx2,zx3,zy1,   &
     &                zy2,zy3,iflag,w,wx,wy,ier)
      implicit none
!*--TVAL4883
      integer iflag , ier
      real x , y , x1 , x2 , x3 , y1 , y2 , y3 , z1 , z2 , z3 , zx1 ,   &
     &     zx2 , zx3 , zy1 , zy2 , zy3 , w , wx , wy
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN FUNCTION VALUES AND FIRST PARTIAL DERIVATIVES AT
! THE THREE VERTICES OF A TRIANGLE, THIS ROUTINE DETERMINES
! A FUNCTION W WHICH AGREES WITH THE GIVEN DATA, RETURNING
! THE VALUE AND (OPTIONALLY) FIRST PARTIAL DERIVATIVES OF W
! AT A POINT (X,Y) IN THE TRIANGLE.  THE INTERPOLATION
! METHOD IS EXACT FOR QUADRATIC POLYNOMIAL DATA.  THE
! TRIANGLE IS PARTITIONED INTO THREE SUBTRIANGLES WITH
! EQUAL AREAS.  W IS CUBIC IN EACH SUBTRIANGLE AND ALONG
! THE EDGES, BUT HAS ONLY ONE CONTINUOUS DERIVATIVE ACROSS
! EDGES.  THE NORMAL DERIVATIVE OF W VARIES LINEARLY ALONG
! EACH OUTER EDGE.  THE VALUES AND PARTIAL DERIVATIVES OF W
! ALONG A TRIANGLE EDGE DEPEND ONLY ON THE DATA VALUES AT
! THE ENDPOINTS OF THE EDGE.  THUS THE METHOD YIELDS C-1
! CONTINUITY WHEN USED TO INTERPOLATE OVER A TRIANGULAR
! GRID.  THIS ALGORITHM IS DUE TO C. L. LAWSON.
!
! INPUT PARAMETERS -   X,Y - COORDINATES OF A POINT AT WHICH
!                            W IS TO BE EVALUATED.
!
!        X1,X2,X3,Y1,Y2,Y3 - COORDINATES OF THE VERTICES OF
!                            A TRIANGLE CONTAINING (X,Y).
!
!                 Z1,Z2,Z3 - FUNCTION VALUES AT THE VERTICES
!                            TO BE INTERPOLATED.
!
!              ZX1,ZX2,ZX3 - X-DERIVATIVE VALUES AT THE
!                            VERTICES.
!
!              ZY1,ZY2,ZY3 - Y-DERIVATIVE VALUES AT THE
!                            VERTICES.
!
!                    IFLAG - OPTION INDICATOR
!                            IFLAG = 0 IF ONLY W IS TO BE
!                                      COMPUTED.
!                            IFLAG = 1 IF W, WX, AND WY ARE
!                                      TO BE RETURNED.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! OUTPUT PARAMETERS -   W - ESTIMATED VALUE OF THE INTERP-
!                           OLATORY FUNCTION AT (X,Y) IF
!                           IER = 0.  OTHERWISE W = 0.
!
!                   WX,WY - PARTIAL DERIVATIVES OF W AT
!                           (X,Y) IF IER = 0 AND IFLAG = 1,
!                           UNCHANGED IF IFLAG .NE. 1, ZERO
!                           IF IER .NE. 0 AND IFLAG = 1.
!
!                     IER - ERROR INDICATOR
!                           IER = 0 IF NO ERRORS WERE
!                                   ENCOUNTERED.
!                           IER = 1 IF THE VERTICES OF THE
!                                   TRIANGLE ARE COLLINEAR.
!
! MODULES REFERENCED BY TVAL - NONE
!
! INTRINSIC FUNCTION CALLED BY TVAL - AMIN1
!
!***********************************************************
!
      integer i , ip1 , ip2 , ip3
      real u(3) , v(3) , sl(3) , area , xp , yp , r(3) , rx(3) , ry(3) ,&
     &     phi(3) , phix(3) , phiy(3) , rmin , c1 , c2 , ro(3) , rox(3) &
     &     , roy(3) , f(3) , g(3) , gx(3) , gy(3) , p(3) , px(3) ,      &
     &     py(3) , q(3) , qx(3) , qy(3) , a(3) , ax(3) , ay(3) , b(3) , &
     &     bx(3) , by(3) , c(3) , cx(3) , cy(3)
!
! LOCAL PARAMETERS -
!
! I =               DO-LOOP INDEX
! IP1,IP2,IP3 =     PERMUTED INDICES FOR COMPUTING RO, ROX,
!                     AND ROY
! U(K) =            X-COMPONENT OF THE VECTOR REPRESENTING
!                     THE SIDE OPPOSITE VERTEX K
! V(K) =            Y-COMPONENT OF THE VECTOR REPRESENTING
!                     THE SIDE OPPOSITE VERTEX K
! SL(K) =           SQUARE OF THE LENGTH OF THE SIDE
!                     OPPOSITE VERTEX K
! AREA =            TWICE THE AREA OF THE TRIANGLE
! XP,YP =           X-X1, Y-Y1
! R(K) =            K-TH BARYCENTRIC COORDINATE
! RX(K),RY(K) =     X,Y PARTIAL DERIVATIVES OF R(K)
! PHI(K)            R(K-1)*R(K+1) -- QUADRATIC
! PHIX(K),PHIY(K) = X,Y PARTIALS OF PHI(K)
! RMIN =            MIN(R1,R2,R3)
! C1,C2 =           FACTORS FOR COMPUTING RO
! RO(K) =           FACTORS FOR COMPUTING G -- CUBIC
!                     CORRECTION TERMS
! ROX(K),ROY(K) =   X,Y PARTIALS OF RO(K)
! F(K) =            FACTORS FOR COMPUTING G, GX, AND GY --
!                     CONSTANT
! G(K) =            FACTORS FOR COMPUTING THE CARDINAL
!                     FUNCTIONS -- CUBIC
! GX(K),GY(K) =     X,Y PARTIALS OF G(K)
! P(K) =            G(K) + PHI(K)
! PX(K),PY(K) =     X,Y PARTIALS OF P(K)
! Q(K) =            G(K) - PHI(K)
! QX(K),QY(K) =     X,Y PARTIALS OF Q(K)
! A(K) =            CARDINAL FUNCTION WHOSE COEFFICIENT IS
!                     Z(K)
! AX(K),AY(K) =     X,Y PARTIALS OF A(K) -- CARDINAL
!                     FUNCTIONS FOR WX AND WY
! B(K) =            TWICE THE CARDINAL FUNCTION WHOSE
!                     COEFFICIENT IS ZX(K)
! BX(K),BY(K) =     X,Y PARTIALS OF B(K)
! C(K) =            TWICE THE CARDINAL FUNCTION WHOSE
!                     COEFFICIENT IS ZY(K)
! CX(K),CY(K) =     X,Y PARTIALS OF C(K)
!
      u(1) = x3 - x2
      u(2) = x1 - x3
      u(3) = x2 - x1
!
      v(1) = y3 - y2
      v(2) = y1 - y3
      v(3) = y2 - y1
!
      do i = 1 , 3
         sl(i) = u(i)*u(i) + v(i)*v(i)
      enddo
!
! AREA = 3-1 X 3-2
!
      area = u(1)*v(2) - u(2)*v(1)
      if ( area==0. ) then
!
! VERTICES ARE COLLINEAR
!
         ier = 1
         w = 0.
         if ( iflag/=1 ) return
         wx = 0.
         wy = 0.
         goto 99999
      else
!
! R(1) = (2-3 X 2-(X,Y))/AREA, R(2) = (1-(X,Y) X 1-3)/AREA,
!   R(3) = (1-2 X 1-(X,Y))/AREA
!
         r(1) = (u(1)*(y-y2)-v(1)*(x-x2))/area
         xp = x - x1
         yp = y - y1
         r(2) = (u(2)*yp-v(2)*xp)/area
         r(3) = (u(3)*yp-v(3)*xp)/area
         ier = 0
!
         phi(1) = r(2)*r(3)
         phi(2) = r(3)*r(1)
         phi(3) = r(1)*r(2)
!
         rmin = amin1(r(1),r(2),r(3))
         if ( rmin==r(1) ) then
            ip1 = 1
            ip2 = 2
            ip3 = 3
         elseif ( rmin/=r(2) ) then
            ip1 = 3
            ip2 = 1
            ip3 = 2
         else
            ip1 = 2
            ip2 = 3
            ip3 = 1
         endif
      endif
!
      c1 = rmin*rmin/2.
      c2 = rmin/3.
      ro(ip1) = (phi(ip1)+5.*c1/3.)*r(ip1) - c1
      ro(ip2) = c1*(r(ip3)-c2)
      ro(ip3) = c1*(r(ip2)-c2)
!
      f(1) = 3.*(sl(2)-sl(3))/sl(1)
      f(2) = 3.*(sl(3)-sl(1))/sl(2)
      f(3) = 3.*(sl(1)-sl(2))/sl(3)
!
      g(1) = (r(2)-r(3))*phi(1) + f(1)*ro(1) - ro(2) + ro(3)
      g(2) = (r(3)-r(1))*phi(2) + f(2)*ro(2) - ro(3) + ro(1)
      g(3) = (r(1)-r(2))*phi(3) + f(3)*ro(3) - ro(1) + ro(2)
!
      do i = 1 , 3
         p(i) = g(i) + phi(i)
         q(i) = g(i) - phi(i)
      enddo
!
      a(1) = r(1) + g(3) - g(2)
      a(2) = r(2) + g(1) - g(3)
      a(3) = r(3) + g(2) - g(1)
!
      b(1) = u(3)*p(3) + u(2)*q(2)
      b(2) = u(1)*p(1) + u(3)*q(3)
      b(3) = u(2)*p(2) + u(1)*q(1)
!
      c(1) = v(3)*p(3) + v(2)*q(2)
      c(2) = v(1)*p(1) + v(3)*q(3)
      c(3) = v(2)*p(2) + v(1)*q(1)
!
! W IS A LINEAR COMBINATION OF THE CARDINAL FUNCTIONS
!
      w = a(1)*z1 + a(2)*z2 + a(3)                                      &
     &    *z3 + (b(1)*zx1+b(2)*zx2+b(3)*zx3+c(1)*zy1+c(2)*zy2+c(3)*zy3) &
     &    /2.
      if ( iflag/=1 ) return
!
! COMPUTE WX AND WY
!
      do i = 1 , 3
         rx(i) = -v(i)/area
         ry(i) = u(i)/area
      enddo
      phix(1) = r(2)*rx(3) + rx(2)*r(3)
      phiy(1) = r(2)*ry(3) + ry(2)*r(3)
      phix(2) = r(3)*rx(1) + rx(3)*r(1)
      phiy(2) = r(3)*ry(1) + ry(3)*r(1)
      phix(3) = r(1)*rx(2) + rx(1)*r(2)
      phiy(3) = r(1)*ry(2) + ry(1)*r(2)
!
      rox(ip1) = rx(ip1)*(phi(ip1)+5.*c1) + r(ip1)*(phix(ip1)-rx(ip1))
      roy(ip1) = ry(ip1)*(phi(ip1)+5.*c1) + r(ip1)*(phiy(ip1)-ry(ip1))
      rox(ip2) = rx(ip1)*(phi(ip2)-c1) + c1*rx(ip3)
      roy(ip2) = ry(ip1)*(phi(ip2)-c1) + c1*ry(ip3)
      rox(ip3) = rx(ip1)*(phi(ip3)-c1) + c1*rx(ip2)
      roy(ip3) = ry(ip1)*(phi(ip3)-c1) + c1*ry(ip2)
!
      gx(1) = (rx(2)-rx(3))*phi(1) + (r(2)-r(3))*phix(1) + f(1)*rox(1)  &
     &        - rox(2) + rox(3)
      gy(1) = (ry(2)-ry(3))*phi(1) + (r(2)-r(3))*phiy(1) + f(1)*roy(1)  &
     &        - roy(2) + roy(3)
      gx(2) = (rx(3)-rx(1))*phi(2) + (r(3)-r(1))*phix(2) + f(2)*rox(2)  &
     &        - rox(3) + rox(1)
      gy(2) = (ry(3)-ry(1))*phi(2) + (r(3)-r(1))*phiy(2) + f(2)*roy(2)  &
     &        - roy(3) + roy(1)
      gx(3) = (rx(1)-rx(2))*phi(3) + (r(1)-r(2))*phix(3) + f(3)*rox(3)  &
     &        - rox(1) + rox(2)
      gy(3) = (ry(1)-ry(2))*phi(3) + (r(1)-r(2))*phiy(3) + f(3)*roy(3)  &
     &        - roy(1) + roy(2)
!
      do i = 1 , 3
         px(i) = gx(i) + phix(i)
         py(i) = gy(i) + phiy(i)
         qx(i) = gx(i) - phix(i)
         qy(i) = gy(i) - phiy(i)
      enddo
!
      ax(1) = rx(1) + gx(3) - gx(2)
      ay(1) = ry(1) + gy(3) - gy(2)
      ax(2) = rx(2) + gx(1) - gx(3)
      ay(2) = ry(2) + gy(1) - gy(3)
      ax(3) = rx(3) + gx(2) - gx(1)
      ay(3) = ry(3) + gy(2) - gy(1)
!
      bx(1) = u(3)*px(3) + u(2)*qx(2)
      by(1) = u(3)*py(3) + u(2)*qy(2)
      bx(2) = u(1)*px(1) + u(3)*qx(3)
      by(2) = u(1)*py(1) + u(3)*qy(3)
      bx(3) = u(2)*px(2) + u(1)*qx(1)
      by(3) = u(2)*py(2) + u(1)*qy(1)
!
      cx(1) = v(3)*px(3) + v(2)*qx(2)
      cy(1) = v(3)*py(3) + v(2)*qy(2)
      cx(2) = v(1)*px(1) + v(3)*qx(3)
      cy(2) = v(1)*py(1) + v(3)*qy(3)
      cx(3) = v(2)*px(2) + v(1)*qx(1)
      cy(3) = v(2)*py(2) + v(1)*qy(1)
!
! WX AND WY ARE LINEAR COMBINATIONS OF THE CARDINAL
!   FUNCTIONS
!
      wx = ax(1)*z1 + ax(2)*z2 + ax(3)                                  &
     &     *z3 + (bx(1)*zx1+bx(2)*zx2+bx(3)*zx3+cx(1)*zy1+cx(2)         &
     &     *zy2+cx(3)*zy3)/2.
      wy = ay(1)*z1 + ay(2)*z2 + ay(3)                                  &
     &     *z3 + (by(1)*zx1+by(2)*zx2+by(3)*zx3+cy(1)*zy1+cy(2)         &
     &     *zy2+cy(3)*zy3)/2.
      return
99999 end subroutine tval
!*==UNIF.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      subroutine unif(n,x,y,z,iadj,iend,nrow,nx,ny,px,py,iflag,zxzy,zz, &
     &                ier)
      implicit none
!*--UNIF5174
      integer n , iadj(*) , iend(n) , nrow , nx , ny , iflag , ier
      real x(n) , y(n) , z(n) , px(nx) , py(ny) , zxzy(2,n) ,           &
     &     zz(nrow,ny)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A THIESSEN TRIANGULATION OF A SET OF POINTS IN THE
! PLANE WITH CORRESPONDING DATA VALUES, THIS ROUTINE INTERP-
! OLATES THE DATA VALUES TO A SET OF RECTANGULAR GRID POINTS
! FOR SUCH APPLICATIONS AS CONTOURING.  THE INTERPOLANT IS
! ONCE CONTINUOUSLY DIFFERENTIABLE.  EXTRAPOLATION IS PER-
! FORMED AT GRID POINTS EXTERIOR TO THE TRIANGULATION.
!
! INPUT PARAMETERS - N - NUMBER OF NODES IN THE TRIANGULA-
!                        TION.  N .GE. 3
!
!                X,Y,Z - N-VECTORS OF NODAL COORDINATES AND
!                        DATA VALUES.
!
!            IADJ,IEND - TRIANGULATION DATA STRUCTURE -- MAY
!                        BE CREATED BY TRMESH.
!
!                 NROW - NUMBER OF ROWS IN THE DIMENSION
!                        STATEMENT OF ZZ.
!
!                NX,NY - NUMBER OF ROWS AND COLUMNS, RESPEC-
!                        TIVELY, IN THE RECTANGULAR GRID.
!                        1 .LE. NX .LE. NROW AND 1 .LE. NY.
!
!                PX,PY - VECTORS OF LENGTH NX AND NY, RE-
!                        SPECTIVELY, CONTAINING THE COORDI-
!                        NATES OF THE GRID LINES.
!
!                IFLAG - OPTION INDICATOR
!                        IFLAG = 0 IF DERIVATIVE ESTIMATES
!                                  AT THE VERTICES OF A
!                                  TRIANGLE ARE TO BE RECOM-
!                                  PUTED FOR EACH GRID POINT
!                                  IN THE TRIANGLE AND NOT
!                                  SAVED.
!                        IFLAG = 1 IF DERIVATIVE ESTIMATES
!                                  ARE INPUT IN ZXZY.
!                        IFLAG = 2 IF DERIVATIVE ESTIMATES
!                                  ARE TO BE COMPUTED ONCE
!                                  FOR EACH NODE (BY GRADL)
!                                  AND SAVED IN ZXZY.
!
!                 ZXZY - NOT USED IF IFLAG = 0, 2 BY N ARRAY
!                        WHOSE COLUMNS CONTAIN THE X AND Y
!                        PARTIAL DERIVATIVE ESTIMATES (X
!                        PARTIALS IN THE FIRST ROW) IF
!                        IFLAG = 1, OR 2 BY N ARRAY (OR WORK
!                        SPACE OF LENGTH .GE. 2*N) IF IFLAG
!                        = 2.
!
! DERIVATIVE ESTIMATES MAY BE COMPUTED BY GRADL OR GRADG.
!
!                   ZZ - NROW BY NCOL ARRAY WITH NROW .GE.
!                        NX AND NCOL .GE. NY.
!
! NONE OF THE INPUT PARAMETERS ARE ALTERED EXCEPT AS
!   NOTED BELOW.
!
! OUTPUT PARAMETERS - ZXZY - 2 BY N ARRAY WHOSE COLUMNS CON-
!                            TAIN X AND Y PARTIAL DERIVATIVE
!                            ESTIMATES AT THE NODES IF IFLAG
!                            = 2 AND IER .GE. 0, NOT ALTERED
!                            IF IFLAG .NE. 2.
!
!                       ZZ - INTERPOLATED VALUES AT THE GRID
!                            POINTS IF IER .GE. 0.
!                            ZZ(I,J) = F(PX(I),PY(J)) FOR
!                            I = 1,...,NX AND J = 1,...,NY.
!
!                      IER - ERROR INDICATOR
!                            IER .GE. 0 IF NO ERRORS WERE
!                                       ENCOUNTERED.  IER
!                                       CONTAINS THE NUMBER
!                                       OF GRID POINTS EXT-
!                                       ERIOR TO THE TRIAN-
!                                       GULATION BOUNDARY.
!                            IER  =  -1 IF N, NX, NY, OR
!                                       IFLAG IS OUT OF
!                                       RANGE.
!                            IER  =  -2 IF THE NODES ARE
!                                       COLLINEAR.
!
! MODULES REFERENCED BY UNIF - INTRC1, TRFIND, TVAL,
!           (AND OPTIONALLY)   GRADL, GETNP, SETUP, GIVENS,
!                                AND ROTATE
!
!***********************************************************
!
      integer nst , ist , nn , ni , nj , ifl , i , j , ierr , nex
      data nst/1/
!
! LOCAL PARAMETERS -
!
! IST =   PARAMETER FOR INTRC1
! NST =   INITIAL VALUE FOR IST
! NN =    LOCAL COPY OF N
! NI,NJ = LOCAL COPIES OF NX AND NY
! IFL =   LOCAL COPY OF IFLAG FOR INTRC1
! I,J =   DO-LOOP INDICES
! IERR =  ERROR FLAG FOR CALLS TO GRADL AND INTRC1
! NEX =   NUMBER OF GRID POINTS EXTERIOR TO THE TRIANGULA-
!           TION BOUNDARY (NUMBER OF EXTRAPOLATED VALUES)
!
      nn = n
      ni = nx
      nj = ny
      ifl = iflag
      if ( nn<3 .or. ni<1 .or. ni>nrow .or. nj<1 .or. ifl<0 .or. ifl>2 )&
     &     then
!
! N, NX, NY, OR IFLAG IS OUT OF RANGE
!
         ier = -1
         return
      else
         ist = nst
         if ( ifl==2 ) then
!
! COMPUTE DERIVATIVE ESTIMATES AT THE NODES.
!
            ifl = 1
            do i = 1 , nn
               call gradl(nn,i,x,y,z,iadj,iend,zxzy(1,i),zxzy(2,i),ierr)
               if ( ierr<0 ) goto 100
            enddo
         endif
!
! COMPUTE INTERPOLATED VALUES
!
         nex = 0
         do j = 1 , nj
            do i = 1 , ni
               call intrc1(nn,px(i),py(j),x,y,z,iadj,iend,ifl,zxzy,ist, &
     &                     zz(i,j),ierr)
               if ( ierr<0 ) goto 100
               nex = nex + ierr
            enddo
         enddo
         ier = nex
         return
      endif
!
! TRIANGULATION NODES ARE COLLINEAR
!
 100  ier = -2
      end subroutine unif
!*==VOLUME.f90  processed by SPAG 7.51DB at 12:04 on  2 Mar 2022
      function volume(n,x,y,z,iadj,iend)
      implicit none
!*--VOLUME5333
!*** Start of declarations inserted by SPAG
      !real trvol , volume
      real volume
!*** End of declarations inserted by SPAG
      integer n , iadj(*) , iend(n)
      real x(n) , y(n) , z(n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       OAK RIDGE NATL. LAB.
!                                             (615) 576-5139
!
!   GIVEN A SET OF N DATA POINTS (X(I),Y(I)) AND FUNCTION
! VALUES Z(I)=F(X(I),Y(I)) AND A TRIANGULATION COVERING THE
! CONVEX HULL H OF THE DATA POINTS, THIS FUNCTION APPROXI-
! MATES THE INTEGRAL OF F OVER H BY INTEGRATING THE PIECE-
! WISE LINEAR INTERPOLANT OF THE DATA VALUES.  VOLUME IS
! PART OF AN INTERPOLATION PACKAGE WHICH PROVIDES ROUTINES
! TO CREATE, UPDATE, AND PLOT THE TRIANGULAR MESH.
!
! INPUT PARAMETERS -      N - NUMBER OF NODES IN THE MESH.
!                             N .GE. 3.
!
!                       X,Y - VECTORS OF COORDINATES OF
!                             THE NODES IN THE MESH.
!
!                         Z - VECTOR OF DATA VALUES AT THE
!                             NODES.
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
! INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.
!
! OUTPUT PARAMETER - VOLUME - SUM OF THE VOLUMES OF THE
!                             LINEAR INTERPOLANTS ON THE
!                             TRIANGLES.
!
! MODULE REFERENCED BY VOLUME - TRVOL
!
!***********************************************************
!
      integer nn , nm2 , n1 , n2 , n3 , indf , indl , indx
      real sum , xn1 , yn1 , zn1
!
! LOCAL PARAMETERS -
!
! NN =          LOCAL COPY OF N
! NM2 =         N-2
! N1,N2,N3 =    VERTICES OF A TRIANGLE IN COUNTERCLOCKWISE
!                 ORDER
! INDF =        IADJ INDEX OF THE FIRST NEIGHBOR OF N1
! INDL =        IADJ INDEX OF THE LAST NEIGHBOR OF N1
! INDX =        IADJ INDEX VARYING FROM INDF TO INDL
! SUM =         TEMPORARY STORAGE FOR ACCUMULATED VOLUME
! XN1,YN1,ZN1 = X(N1), Y(N1), Z(N1) -- STORED LOCALLY FOR
!                 EFFICIENCY
!
      nn = n
      if ( nn<3 ) then
!
! N IS OUT OF RANGE
!
         volume = 0.
         goto 99999
      endif
!
! INITIALIZATION
!
      nm2 = nn - 2
      indf = 1
      sum = 0.
!
! LOOP ON TRIANGLES (N1,N2,N3) SUCH THAT N2 AND N3 ARE
!   ADJACENT NEIGHBORS OF N1 WHICH ARE BOTH LARGER THAN N1
!
      do n1 = 1 , nm2
         xn1 = x(n1)
         yn1 = y(n1)
         zn1 = z(n1)
         indl = iend(n1)
         do indx = indf , indl
            n2 = iadj(indx)
            n3 = iadj(indx+1)
            if ( indx==indl ) n3 = iadj(indf)
            if ( n2>=n1 .and. n3>=n1 ) sum = sum +                      &
     &           trvol(xn1,x(n2),x(n3),yn1,y(n2),y(n3),zn1,z(n2),z(n3))
         enddo
         indf = indl + 1
      enddo
!
      volume = sum
      return
99999 end function volume

      subroutine circum (x1,x2,x3,y1,y2,y3, cx,cy,ier)
      integer ier
      real    x1, x2, x3, y1, y2, y3, cx, cy
!
!***********************************************************
!
!                                               robert renka
!                                       oak ridge natl. lab.
!                                             (615) 576-5139
!
!   this subroutine computes the coordinates of the center
! of a circle defined by three points in the plane.
!
! input parameters -
!
!       x1,...,y3 - x and y coordinates of three points in
!                   the plane.
!
! input parameters are not altered by this routine.
!
! output parameters -
!
!       cx,cy - coordinates of the center of the circle
!               unless ier = 1.
!
!       ier - error indicator
!             ier = 0 if no errors were encountered.
!             ier = 1 if the points are collinear.
!
! modules required by circum - none
!
!***********************************************************
!
      real u(3), v(3), ds(3)
!
! set u(k) and v(k) to the x and y components of the edge
!   opposite vertex k, treating the points as vertices of
!   a triangle.
!
      u(1) = x3 - x2
      u(2) = x1 - x3
      u(3) = x2 - x1
      v(1) = y3 - y2
      v(2) = y1 - y3
      v(3) = y2 - y1
!
! set a to twice the signed area of the triangle.  a .gt. 0
!   iff (x3,y3) is strictly to the left of the edge from
!   (x1,y1) to (x2,y2).
!
      a = u(1)*v(2) - u(2)*v(1)

      if ( a == 0.0 ) then
!
! collinear points
!
         ier = 1
         return
      endif
!
! set ds(k) to the squared distance from the origin to
!   vertex k.
!
      ds(1) = x1**2 + y1**2
      ds(2) = x2**2 + y2**2
      ds(3) = x3**2 + y3**2
!
! compute factors of cx and cy.
!
      fx = 0.
      fy = 0.
      do i = 1,3
        fx = fx - ds(i)*v(i)
        fy = fy + ds(i)*u(i)
      enddo
      cx = fx/2./a
      cy = fy/2./a
      ier = 0
      end subroutine circum

end module alg624
