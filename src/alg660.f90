module alg660
    use ieee_arithmetic
contains
!*==QSHEP2.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine qshep2(n,x,y,f,nq,nw,nr,lcell,lnext,xmin,ymin,dx,dy,   &
     &                  rmax,rsq,a,ier)
      implicit none
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - F
!A INPUT  - NQ
!A INPUT  - NW
!A INPUT  - NR
!A PASSED - LCELL
!A OUTPUT - LNEXT
!A OUTPUT - XMIN
!A OUTPUT - YMIN
!A OUTPUT - DX
!A OUTPUT - DY
!A OUTPUT - RMAX
!A OUTPUT - RSQ
!A OUTPUT - A
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       GETNP2   GIVENS   ROTATE   SETUP2   STORE2
! called by   AA0005
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  AV       AVSQ     B        C        DDX      DDY      DMIN     DTOL     FK       I        IB       IERR     IP1      IRM1     IROW     J        JP1      K        LMAX     LNP      NEQ      NN       NNQ      NNR      NNW      NP       NPTS
!             NQWMAX   RQ       RS       RSMX     RSOLD    RTOL     RWS      S        SF       SUM      T        XK       XMN      YK       YMN
! uses PARAMs *** NONE ****
!*++********************************************************************
      integer n , nq , nw , nr , lcell(nr,nr) , lnext(n) , ier
      real x(n) , y(n) , f(n) , xmin , ymin , dx , dy , rmax , rsq(n) , &
     &     a(5,n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!                                                   01/08/90
!
!   THIS SUBROUTINE COMPUTES A SET OF PARAMETERS A AND RSQ
! DEFINING A SMOOTH (ONCE CONTINUOUSLY DIFFERENTIABLE) BI-
! VARIATE FUNCTION Q(X,Y) WHICH INTERPOLATES DATA VALUES F
! AT SCATTERED NODES (X,Y).  THE INTERPOLANT Q MAY BE EVAL-
! UATED AT AN ARBITRARY POINT BY FUNCTION QS2VAL, AND ITS
! FIRST DERIVATIVES ARE COMPUTED BY SUBROUTINE QS2GRD.
!   THE INTERPOLATION SCHEME IS A MODIFIED QUADRATIC SHEPARD
! METHOD --
!
! Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))
!
! FOR BIVARIATE FUNCTIONS W(K) AND Q(K).  THE NODAL FUNC-
! TIONS ARE GIVEN BY
!
! Q(K)(X,Y) = A(1,K)*(X-X(K))**2 + A(2,K)*(X-X(K))*(Y-Y(K))
!           + A(3,K)*(Y-Y(K))**2 + A(4,K)*(X-X(K))
!           + A(5,K)*(Y-Y(K))    + F(K) .
!
! THUS, Q(K) IS A QUADRATIC FUNCTION WHICH INTERPOLATES THE
! DATA VALUE AT NODE K.  ITS COEFFICIENTS A(,K) ARE OBTAINED
! BY A WEIGHTED LEAST SQUARES FIT TO THE CLOSEST NQ DATA
! POINTS WITH WEIGHTS SIMILAR TO W(K).  NOTE THAT THE RADIUS
! OF INFLUENCE FOR THE LEAST SQUARES FIT IS FIXED FOR EACH
! K, BUT VARIES WITH K.
!   THE WEIGHTS ARE TAKEN TO BE
!
! W(K)(X,Y) = ( (R(K)-D(K))+ / R(K)*D(K) )**2
!
! WHERE (R(K)-D(K))+ = 0 IF R(K) .LE. D(K) AND D(K)(X,Y) IS
! THE EUCLIDEAN DISTANCE BETWEEN (X,Y) AND (X(K),Y(K)).  THE
! RADIUS OF INFLUENCE R(K) VARIES WITH K AND IS CHOSEN SO
! THAT NW NODES ARE WITHIN THE RADIUS.  NOTE THAT W(K) IS
! NOT DEFINED AT NODE (X(K),Y(K)), BUT Q(X,Y) HAS LIMIT F(K)
! AS (X,Y) APPROACHES (X(K),Y(K)).
!
! ON INPUT --
!
!       N = NUMBER OF NODES AND ASSOCIATED DATA VALUES.
!           N .GE. 6.
!
!       X,Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!             COORDINATES OF THE NODES.
!
!       F = ARRAY OF LENGTH N CONTAINING THE DATA VALUES
!           IN ONE-TO-ONE CORRESPONDENCE WITH THE NODES.
!
!       NQ = NUMBER OF DATA POINTS TO BE USED IN THE LEAST
!            SQUARES FIT FOR COEFFICIENTS DEFINING THE NODAL
!            FUNCTIONS Q(K).  A HIGHLY RECOMMENDED VALUE IS
!            NQ = 13.  5 .LE. NQ .LE. MIN(40,N-1).
!
!       NW = NUMBER OF NODES WITHIN (AND DEFINING) THE RADII
!            OF INFLUENCE R(K) WHICH ENTER INTO THE WEIGHTS
!            W(K).  FOR N SUFFICIENTLY LARGE, A RECOMMENDED
!            VALUE IS NW = 19.  1 .LE. NW .LE. MIN(40,N-1).
!
!       NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID DE-
!            FINED IN SUBROUTINE STORE2.  A RECTANGLE CON-
!            TAINING THE NODES IS PARTITIONED INTO CELLS IN
!            ORDER TO INCREASE SEARCH EFFICIENCY.  NR =
!            SQRT(N/3) IS RECOMMENDED.  NR .GE. 1.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!       LCELL = ARRAY OF LENGTH .GE. NR**2.
!
!       LNEXT = ARRAY OF LENGTH .GE. N.
!
!       RSQ = ARRAY OF LENGTH .GE. N.
!
!       A = ARRAY OF LENGTH .GE. 5N.
!
! ON OUTPUT --
!
!       LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED
!               WITH CELLS.  REFER TO STORE2.
!
!       LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!               CES.  REFER TO STORE2.
!
!       XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL
!                         DIMENSIONS.  REFER TO STORE2.
!
!       RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!              MAXIMUM RADIUS R(K).
!
!       RSQ = ARRAY CONTAINING THE SQUARES OF THE RADII R(K)
!             WHICH ENTER INTO THE WEIGHTS W(K).
!
!       A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR
!           QUADRATIC NODAL FUNCTION Q(K) IN COLUMN K.
!
!   NOTE THAT THE ABOVE OUTPUT PARAMETERS ARE NOT DEFINED
! UNLESS IER = 0.
!
!       IER = ERROR INDICATOR --
!             IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!             IER = 1 IF N, NQ, NW, OR NR IS OUT OF RANGE.
!             IER = 2 IF DUPLICATE NODES WERE ENCOUNTERED.
!             IER = 3 IF ALL NODES ARE COLLINEAR.
!
! MODULES REQUIRED BY QSHEP2 -- GETNP2, GIVENS, ROTATE,
!                                 SETUP2, STORE2
!
! INTRINSIC FUNCTIONS CALLED BY QSHEP2 -- ABS, AMIN1, FLOAT,
!                                           MAX0, MIN0, SQRT
!
!***********************************************************
!
      integer i , ib , ierr , ip1 , irm1 , irow , j , jp1 , k , lmax ,  &
     &        lnp , neq , nn , nnq , nnr , nnw , np , npts(40) , nqwmax
      real av , avsq , b(6,6) , c , ddx , ddy , dmin , dtol , fk , rq , &
     &     rs , rsmx , rsold , rtol , rws , s , sf , sum , t , xk ,     &
     &     xmn , yk , ymn
!
      data rtol/1.E-5/ , dtol/.01/ , sf/1./
!
! LOCAL PARAMETERS --
!
! AV =         ROOT-MEAN-SQUARE DISTANCE BETWEEN K AND THE
!                NODES IN THE LEAST SQUARES SYSTEM (UNLESS
!                ADDITIONAL NODES ARE INTRODUCED FOR STABIL-
!                ITY).  THE FIRST 3 COLUMNS OF THE MATRIX
!                ARE SCALED BY 1/AVSQ, THE LAST 2 BY 1/AV
! AVSQ =       AV*AV
! B =          TRANSPOSE OF THE AUGMENTED REGRESSION MATRIX
! C =          FIRST COMPONENT OF THE PLANE ROTATION USED TO
!                ZERO THE LOWER TRIANGLE OF B**T -- COMPUTED
!                BY SUBROUTINE GIVENS
! DDX,DDY =    LOCAL VARIABLES FOR DX AND DY
! DMIN =       MINIMUM OF THE MAGNITUDES OF THE DIAGONAL
!                ELEMENTS OF THE REGRESSION MATRIX AFTER
!                ZEROS ARE INTRODUCED BELOW THE DIAGONAL
! DTOL =       TOLERANCE FOR DETECTING AN ILL-CONDITIONED
!                SYSTEM.  THE SYSTEM IS ACCEPTED WHEN DMIN
!                .GE. DTOL
! FK =         DATA VALUE AT NODE K -- F(K)
! I =          INDEX FOR A, B, AND NPTS
! IB =         DO-LOOP INDEX FOR BACK SOLVE
! IERR =       ERROR FLAG FOR THE CALL TO STORE2
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       ROW INDEX FOR B
! J =          INDEX FOR A AND B
! JP1 =        J+1
! K =          NODAL FUNCTION INDEX AND COLUMN INDEX FOR A
! LMAX =       MAXIMUM NUMBER OF NPTS ELEMENTS (MUST BE CON-
!                SISTENT WITH THE DIMENSION STATEMENT ABOVE)
! LNP =        CURRENT LENGTH OF NPTS
! NEQ =        NUMBER OF EQUATIONS IN THE LEAST SQUARES FIT
! NN,NNQ,NNR = LOCAL COPIES OF N, NQ, AND NR
! NNW =        LOCAL COPY OF NW
! NP =         NPTS ELEMENT
! NPTS =       ARRAY CONTAINING THE INDICES OF A SEQUENCE OF
!                NODES TO BE USED IN THE LEAST SQUARES FIT
!                OR TO COMPUTE RSQ.  THE NODES ARE ORDERED
!                BY DISTANCE FROM K AND THE LAST ELEMENT
!                (USUALLY INDEXED BY LNP) IS USED ONLY TO
!                DETERMINE RQ, OR RSQ(K) IF NW .GT. NQ
! NQWMAX =     MAX(NQ,NW)
! RQ =         RADIUS OF INFLUENCE WHICH ENTERS INTO THE
!                WEIGHTS FOR Q(K) (SEE SUBROUTINE SETUP2)
! RS =         SQUARED DISTANCE BETWEEN K AND NPTS(LNP) --
!                USED TO COMPUTE RQ AND RSQ(K)
! RSMX =       MAXIMUM RSQ ELEMENT ENCOUNTERED
! RSOLD =      SQUARED DISTANCE BETWEEN K AND NPTS(LNP-1) --
!                USED TO COMPUTE A RELATIVE CHANGE IN RS
!                BETWEEN SUCCEEDING NPTS ELEMENTS
! RTOL =       TOLERANCE FOR DETECTING A SUFFICIENTLY LARGE
!                RELATIVE CHANGE IN RS.  IF THE CHANGE IS
!                NOT GREATER THAN RTOL, THE NODES ARE
!                TREATED AS BEING THE SAME DISTANCE FROM K
! RWS =        CURRENT VALUE OF RSQ(K)
! S =          SECOND COMPONENT OF THE PLANE GIVENS ROTATION
! SF =         MARQUARDT STABILIZATION FACTOR USED TO DAMP
!                OUT THE FIRST 3 SOLUTION COMPONENTS (SECOND
!                PARTIALS OF THE QUADRATIC) WHEN THE SYSTEM
!                IS ILL-CONDITIONED.  AS SF INCREASES, THE
!                FITTING FUNCTION APPROACHES A LINEAR
! SUM =        SUM OF SQUARED EUCLIDEAN DISTANCES BETWEEN
!                NODE K AND THE NODES USED IN THE LEAST
!                SQUARES FIT (UNLESS ADDITIONAL NODES ARE
!                ADDED FOR STABILITY)
! T =          TEMPORARY VARIABLE FOR ACCUMULATING A SCALAR
!                PRODUCT IN THE BACK SOLVE
! XK,YK =      COORDINATES OF NODE K -- X(K), Y(K)
! XMN,YMN =    LOCAL VARIABLES FOR XMIN AND YMIN
!
      nn = n
      nnq = nq
      nnw = nw
      nnr = nr
      nqwmax = max0(nnq,nnw)
      lmax = min0(40,nn-1)
      if ( 5>nnq .or. 1>nnw .or. nqwmax>lmax .or. nnr<1 ) then
!
! N, NQ, NW, OR NR IS OUT OF RANGE.
!
         ier = 1
         return
      else
!
! CREATE THE CELL DATA STRUCTURE, AND INITIALIZE RSMX.
!
         call store2(nn,x,y,nnr,lcell,lnext,xmn,ymn,ddx,ddy,ierr)
         if ( ierr/=0 ) goto 100
         rsmx = 0.
!
! OUTER LOOP ON NODE K
!
         do k = 1 , nn
            xk = x(k)
            yk = y(k)
            fk = f(k)
!
! MARK NODE K TO EXCLUDE IT FROM THE SEARCH FOR NEAREST
!   NEIGHBORS.
!
            lnext(k) = -lnext(k)
!
! INITIALIZE FOR LOOP ON NPTS.
!
            rs = 0.
            sum = 0.
            rws = 0.
            rq = 0.
            lnp = 0
            do
!
! COMPUTE NPTS, LNP, RWS, NEQ, RQ, AND AVSQ.
!
               sum = sum + rs
               if ( lnp==lmax ) then
!
! ALL LMAX NODES ARE INCLUDED IN NPTS.  RWS AND/OR RQ**2 IS
!   (ARBITRARILY) TAKEN TO BE 10 PERCENT LARGER THAN THE
!   DISTANCE RS TO THE LAST NODE INCLUDED.
!
                  if ( rws==0. ) rws = 1.1*rs
                  if ( rq==0. ) then
                     neq = lmax
                     rq = sqrt(1.1*rs)
                     avsq = sum/float(neq)
                  endif
                  exit
               else
                  lnp = lnp + 1
                  rsold = rs
                  call getnp2(xk,yk,x,y,nnr,lcell,lnext,xmn,ymn,ddx,ddy,&
     &                        np,rs)
                  if ( rs==0. ) goto 50
                  npts(lnp) = np
                  if ( (rs-rsold)/rs>=rtol ) then
                     if ( rws==0. .and. lnp>nnw ) rws = rs
                     if ( rq==0. .and. lnp>nnq ) then
!
!   RQ = 0 (NOT YET COMPUTED) AND LNP .GT. NQ.  RQ =
!     SQRT(RS) IS SUFFICIENTLY LARGE TO (STRICTLY) INCLUDE
!     NQ NODES.  THE LEAST SQUARES FIT WILL INCLUDE NEQ =
!     LNP - 1 EQUATIONS FOR 5 .LE. NQ .LE. NEQ .LT. LMAX
!     .LE. N-1.
!
                        neq = lnp - 1
                        rq = sqrt(rs)
                        avsq = sum/float(neq)
                     endif
!
!   BOTTOM OF LOOP -- TEST FOR TERMINATION.
!
                     if ( lnp>nqwmax ) exit
                  endif
               endif
            enddo
!
! STORE RSQ(K), UPDATE RSMX IF NECESSARY, AND COMPUTE AV.
!
            rsq(k) = rws
            if ( rws>rsmx ) rsmx = rws
            av = sqrt(avsq)
!
! SET UP THE AUGMENTED REGRESSION MATRIX (TRANSPOSED) AS THE
!   COLUMNS OF B, AND ZERO OUT THE LOWER TRIANGLE (UPPER
!   TRIANGLE OF B) WITH GIVENS ROTATIONS -- QR DECOMPOSITION
!   WITH ORTHOGONAL MATRIX Q NOT STORED.
!
            i = 0
            do
               i = i + 1
               np = npts(i)
               irow = min0(i,6)
               call setup2(xk,yk,fk,x(np),y(np),f(np),av,avsq,rq,       &
     &                     b(1,irow))
               if ( i/=1 ) then
                  irm1 = irow - 1
                  do j = 1 , irm1
                     jp1 = j + 1
                     call givens(b(j,j),b(j,irow),c,s)
                     call rotate(6-j,c,s,b(jp1,j),b(jp1,irow))
                  enddo
                  if ( i>=neq ) then
!
! TEST THE SYSTEM FOR ILL-CONDITIONING.
!
                     dmin = amin1(abs(b(1,1)),abs(b(2,2)),abs(b(3,3)),  &
     &                      abs(b(4,4)),abs(b(5,5)))
                     if ( dmin*rq<dtol ) then
                        if ( neq==lmax ) then
!
! STABILIZE THE SYSTEM BY DAMPING SECOND PARTIALS -- ADD
!   MULTIPLES OF THE FIRST THREE UNIT VECTORS TO THE FIRST
!   THREE EQUATIONS.
!
                           do i = 1 , 3
                              b(i,6) = sf
                              ip1 = i + 1
                              do j = ip1 , 6
                                 b(j,6) = 0.
                              enddo
                              do j = i , 5
                                 jp1 = j + 1
                                 call givens(b(j,j),b(j,6),c,s)
                                 call rotate(6-j,c,s,b(jp1,j),b(jp1,6))
                              enddo
                           enddo
!
! TEST THE STABILIZED SYSTEM FOR ILL-CONDITIONING.
!
                           dmin = amin1(abs(b(1,1)),abs(b(2,2)),        &
     &                            abs(b(3,3)),abs(b(4,4)),abs(b(5,5)))
                           if ( dmin*rq<dtol ) goto 100
                        else
                           do
!
! INCREASE RQ AND ADD ANOTHER EQUATION TO THE SYSTEM TO
!   IMPROVE THE CONDITIONING.  THE NUMBER OF NPTS ELEMENTS
!   IS ALSO INCREASED IF NECESSARY.
!
                              rsold = rs
                              neq = neq + 1
                              if ( neq==lmax ) then
!
                                 rq = sqrt(1.1*rs)
                                 goto 20
                              elseif ( neq==lnp ) then
!
!   ADD AN ELEMENT TO NPTS.
!
                                 lnp = lnp + 1
                                 call getnp2(xk,yk,x,y,nnr,lcell,lnext, &
     &                              xmn,ymn,ddx,ddy,np,rs)
                                 if ( np==0 ) goto 50
                                 npts(lnp) = np
                                 if ( (rs-rsold)/rs>=rtol ) then
                                    rq = sqrt(rs)
                                    goto 20
                                 endif
                              else
!
!   NEQ .LT. LNP
!
                                 np = npts(neq+1)
                                 rs = (x(np)-xk)**2 + (y(np)-yk)**2
                                 if ( (rs-rsold)/rs>=rtol ) then
                                    rq = sqrt(rs)
                                    goto 20
                                 endif
                              endif
                           enddo
                        endif
                     endif
!
! SOLVE THE 5 BY 5 TRIANGULAR SYSTEM FOR THE COEFFICIENTS
!
                     do ib = 1 , 5
                        i = 6 - ib
                        t = 0.
                        if ( i/=5 ) then
                           ip1 = i + 1
                           do j = ip1 , 5
                              t = t + b(j,i)*a(j,k)
                           enddo
                        endif
                        a(i,k) = (b(6,i)-t)/b(i,i)
                     enddo
!
! SCALE THE COEFFICIENTS TO ADJUST FOR THE COLUMN SCALING.
!
                     do i = 1 , 3
                        a(i,k) = a(i,k)/avsq
                     enddo
                     a(4,k) = a(4,k)/av
                     a(5,k) = a(5,k)/av
!
! UNMARK K AND THE ELEMENTS OF NPTS.
!
                     lnext(k) = -lnext(k)
                     do i = 1 , lnp
                        np = npts(i)
                        lnext(np) = -lnext(np)
                     enddo
                     exit
                  endif
               endif
 20         enddo
         enddo
!
! NO ERRORS ENCOUNTERED.
!
         xmin = xmn
         ymin = ymn
         dx = ddx
         dy = ddy
         rmax = sqrt(rsmx)
         ier = 0
         return
!
! DUPLICATE NODES WERE ENCOUNTERED BY GETNP2.
!
 50      ier = 2
         return
      endif
!
! NO UNIQUE SOLUTION DUE TO COLLINEAR NODES.
!
 100  xmin = xmn
      ymin = ymn
      dx = ddx
      dy = ddy
      ier = 3
      end subroutine qshep2
!*==QS2VAL.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      function qs2val(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,rmax,&
     &                rsq,a,ier)
      implicit none
!*--********************************************************************
!A INPUT  - PX
!A INPUT  - PY
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - F
!A INPUT  - NR
!A INPUT  - LCELL
!A INPUT  - LNEXT
!A INPUT  - XMIN
!A INPUT  - YMIN
!A INPUT  - DX
!A INPUT  - DY
!A INPUT  - RMAX
!A INPUT  - RSQ
!A INPUT  - A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0005
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DELX     DELY     DS       DXSQ     DYSQ     I        IMAX     IMIN     J        JMAX     JMIN     K        KP       RD       RDS      RS       SW       SWQ      W        XP       YP
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real delx , dely , ds , dxsq , dysq , qs2val , rd , rds , rs ,    &
     &     sw , swq , w , xp , yp
      integer i , imax , imin , j , jmax , jmin , k , kp , ier
!*** End of declarations inserted by SPAG
      integer n , nr , lcell(nr,nr) , lnext(n)
      real px , py , x(n) , y(n) , f(n) , xmin , ymin , dx , dy , rmax ,&
     &     rsq(n) , a(5,n)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!                                                   10/28/87
!
!   THIS FUNCTION RETURNS THE VALUE Q(PX,PY) WHERE Q IS THE
! WEIGHTED SUM OF QUADRATIC NODAL FUNCTIONS DEFINED IN SUB-
! ROUTINE QSHEP2.  QS2GRD MAY BE CALLED TO COMPUTE A GRADI-
! ENT OF Q ALONG WITH THE VALUE, AND/OR TO TEST FOR ERRORS.
!
! ON INPUT --
!
!       PX,PY = CARTESIAN COORDINATES OF THE POINT P AT
!               WHICH Q IS TO BE EVALUATED.
!
!       N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!           N .GE. 6.
!
!       X,Y,F = ARRAYS OF LENGTH N CONTAINING THE NODES AND
!               DATA VALUES INTERPOLATED BY Q.
!
!       NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!            REFER TO STORE2.  NR .GE. 1.
!
!       LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED
!               WITH CELLS.  REFER TO STORE2.
!
!       LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!               CES.  REFER TO STORE2.
!
!       XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL
!                         DIMENSIONS.  DX AND DY MUST BE
!                         POSITIVE.  REFER TO STORE2.
!
!       RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!              MAXIMUM RADIUS.
!
!       RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!             WHICH ENTER INTO THE WEIGHTS DEFINING Q.
!
!       A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!           NODAL FUNCTIONS DEFINING Q.
!
!   INPUT PARAMETERS ARE NOT ALTERED BY THIS FUNCTION.  THE
! PARAMETERS OTHER THAN PX AND PY SHOULD BE INPUT UNALTERED
! FROM THEIR VALUES ON OUTPUT FROM QSHEP2.  THIS FUNCTION
! SHOULD NOT BE CALLED IF A NONZERO ERROR FLAG WAS RETURNED
! BY QSHEP2.
!
! ON OUTPUT --
!
!       QS2VAL = FUNCTION VALUE Q(PX,PY) UNLESS N, NR, DX,
!                DY, OR RMAX IS INVALID, IN WHICH CASE NO
!                VALUE IS RETURNED.
!
! MODULES REQUIRED BY QS2VAL -- NONE
!
! INTRINSIC FUNCTIONS CALLED BY QS2VAL -- IFIX, SQRT
!
!***********************************************************
!
      ier = 0
      xp  = px
      yp  = py
      if ( n<6 .or. nr<1 .or. dx<=0. .or. dy<=0. .or. rmax<0. )then
          ier    = 1
          qs2val = ieee_value( qs2val, ieee_quiet_nan )
          return
      endif
!
! SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING
!   THE RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE
!   P.  THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTER-
!   SECTED BY (OR CONTAINED IN) A CIRCLE OF RADIUS RMAX
!   CENTERED AT P.
!
      imin = ifix((xp-xmin-rmax)/dx) + 1
      imax = ifix((xp-xmin+rmax)/dx) + 1
      if ( imin<1 ) imin = 1
      if ( imax>nr ) imax = nr
      jmin = ifix((yp-ymin-rmax)/dy) + 1
      jmax = ifix((yp-ymin+rmax)/dy) + 1
      if ( jmin<1 ) jmin = 1
      if ( jmax>nr ) jmax = nr
!
! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE
!   OF RADIUS RMAX.
!
      if ( imin>imax .or. jmin>jmax ) then
!
! ALL WEIGHTS ARE 0 AT P.
!
         ier    = 2
         qs2val = ieee_value( qs2val, ieee_quiet_nan )
      else
!
! ACCUMULATE WEIGHT VALUES IN SW AND WEIGHTED NODAL FUNCTION
!   VALUES IN SWQ.  THE WEIGHTS ARE W(K) = ((R-D)+/(R*D))**2
!   FOR R**2 = RSQ(K) AND D = DISTANCE BETWEEN P AND NODE K.
!
         sw = 0.
         swq = 0.
!
! OUTER LOOP ON CELLS (I,J).
!
         do j = jmin , jmax
            do i = imin , imax
               k = lcell(i,j)
               if ( k==0 ) cycle
!
! INNER LOOP ON NODES K.
!
 10            delx = xp - x(k)
               dely = yp - y(k)
               dxsq = delx*delx
               dysq = dely*dely
               ds = dxsq + dysq
               rs = rsq(k)
               if ( ds<rs ) then
                  if ( ds==0. ) goto 50
                  rds = rs*ds
                  rd = sqrt(rds)
                  w = (rs+ds-rd-rd)/rds
                  sw = sw + w
                  swq = swq + w*(a(1,k)*dxsq+a(2,k)*delx*dely+a(3,k)    &
     &                  *dysq+a(4,k)*delx+a(5,k)*dely+f(k))
               endif
!
! BOTTOM OF LOOP ON NODES IN CELL (I,J).
!
               kp = k
               k = lnext(kp)
               if ( k/=kp ) goto 10
            enddo
         enddo
!
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.
!
         if ( sw==0. ) then
            ier    = 2
            qs2val = ieee_value( qs2val, ieee_quiet_nan )
            goto 99999
         else
            qs2val = swq/sw
            return
         endif
!
! (PX,PY) = (X(K),Y(K))
!
 50      qs2val = f(k)
         return
      endif
99999 end function qs2val
!*==QS2GRD.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine qs2grd(px,py,n,x,y,f,nr,lcell,lnext,xmin,ymin,dx,dy,   &
     &                  rmax,rsq,a,q,qx,qy,ier)
      implicit none
!*--********************************************************************
!A INPUT  - PX
!A INPUT  - PY
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - F
!A INPUT  - NR
!A INPUT  - LCELL
!A INPUT  - LNEXT
!A INPUT  - XMIN
!A INPUT  - YMIN
!A INPUT  - DX
!A INPUT  - DY
!A INPUT  - RMAX
!A INPUT  - RSQ
!A INPUT  - A
!A OUTPUT - Q
!A OUTPUT - QX
!A OUTPUT - QY
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   AA0005
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DELX     DELY     DS       DXSQ     DYSQ     I        IMAX     IMIN     J        JMAX     JMIN     K        KP       QK       QKX      QKY      RD       RDS      RS       SW       SWQ      SWQX     SWQY     SWS      SWX      SWY      T
!             W        WX       WY       XP       YP
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real delx , dely , ds , dxsq , dysq , qk , qkx , qky , rd , rds , &
     &     rs , sw , swq , swqx , swqy , sws , swx , swy , t , w
      real wx , wy , xp , yp
      integer i , imax , imin , j , jmax , jmin , k , kp
!*** End of declarations inserted by SPAG
      integer n , nr , lcell(nr,nr) , lnext(n) , ier
      real px , py , x(n) , y(n) , f(n) , xmin , ymin , dx , dy , rmax ,&
     &     rsq(n) , a(5,n) , q , qx , qy
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!                                                   10/28/87
!
!   THIS SUBROUTINE COMPUTES THE VALUE AND GRADIENT AT
! (PX,PY) OF THE INTERPOLATORY FUNCTION Q DEFINED IN SUB-
! ROUTINE QSHEP2.  Q(X,Y) IS A WEIGHTED SUM OF QUADRATIC
! NODAL FUNCTIONS.
!
! ON INPUT --
!
!       PX,PY = CARTESIAN COORDINATES OF THE POINT AT WHICH
!               Q AND ITS PARTIALS ARE TO BE EVALUATED.
!
!       N = NUMBER OF NODES AND DATA VALUES DEFINING Q.
!           N .GE. 6.
!
!       X,Y,F = ARRAYS OF LENGTH N CONTAINING THE NODES AND
!               DATA VALUES INTERPOLATED BY Q.
!
!       NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!            REFER TO STORE2.  NR .GE. 1.
!
!       LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED
!               WITH CELLS.  REFER TO STORE2.
!
!       LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!               CES.  REFER TO STORE2.
!
!       XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL
!                         DIMENSIONS.  DX AND DY MUST BE
!                         POSITIVE.  REFER TO STORE2.
!
!       RMAX = SQUARE ROOT OF THE LARGEST ELEMENT IN RSQ --
!              MAXIMUM RADIUS.
!
!       RSQ = ARRAY OF LENGTH N CONTAINING THE SQUARED RADII
!             WHICH ENTER INTO THE WEIGHTS DEFINING Q.
!
!       A = 5 BY N ARRAY CONTAINING THE COEFFICIENTS FOR THE
!           NODAL FUNCTIONS DEFINING Q.
!
!   INPUT PARAMETERS ARE NOT ALTERED BY THIS SUBROUTINE.
! THE PARAMETERS OTHER THAN PX AND PY SHOULD BE INPUT UNAL-
! TERED FROM THEIR VALUES ON OUTPUT FROM QSHEP2.  THIS SUB-
! ROUTINE SHOULD NOT BE CALLED IF A NONZERO ERROR FLAG WAS
! RETURNED BY QSHEP2.
!
! ON OUTPUT --
!
!       Q = VALUE OF Q AT (PX,PY) UNLESS IER .EQ. 1, IN
!           WHICH CASE NO VALUES ARE RETURNED.
!
!       QX,QY = FIRST PARTIAL DERIVATIVES OF Q AT (PX,PY)
!               UNLESS IER .EQ. 1.
!
!       IER = ERROR INDICATOR
!             IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!             IER = 1 IF N, NR, DX, DY OR RMAX IS INVALID.
!             IER = 2 IF NO ERRORS WERE ENCOUNTERED BUT
!                     (PX,PY) IS NOT WITHIN THE RADIUS R(K)
!                     FOR ANY NODE K (AND THUS Q=QX=QY=0).
!
! MODULES REQUIRED BY QS2GRD -- NONE
!
! INTRINSIC FUNCTIONS CALLED BY QS2GRD -- IFIX, SQRT
!
!***********************************************************
!
      xp = px
      yp = py
      if ( n<6 .or. nr<1 .or. dx<=0. .or. dy<=0. .or. rmax<0. ) then
!
! INVALID INPUT PARAMETER.
!
         ier = 1
         return
      else
!
! SET IMIN, IMAX, JMIN, AND JMAX TO CELL INDICES DEFINING
!   THE RANGE OF THE SEARCH FOR NODES WHOSE RADII INCLUDE
!   P.  THE CELLS WHICH MUST BE SEARCHED ARE THOSE INTER-
!   SECTED BY (OR CONTAINED IN) A CIRCLE OF RADIUS RMAX
!   CENTERED AT P.
!
         imin = ifix((xp-xmin-rmax)/dx) + 1
         imax = ifix((xp-xmin+rmax)/dx) + 1
         if ( imin<1 ) imin = 1
         if ( imax>nr ) imax = nr
         jmin = ifix((yp-ymin-rmax)/dy) + 1
         jmax = ifix((yp-ymin+rmax)/dy) + 1
         if ( jmin<1 ) jmin = 1
         if ( jmax>nr ) jmax = nr
!
! THE FOLLOWING IS A TEST FOR NO CELLS WITHIN THE CIRCLE
!   OF RADIUS RMAX.
!
         if ( imin<=imax .and. jmin<=jmax ) then
!
! Q = SWQ/SW = SUM(W(K)*Q(K))/SUM(W(K)) WHERE THE SUM IS
!   FROM K = 1 TO N, Q(K) IS THE QUADRATIC NODAL FUNCTION,
!   AND W(K) = ((R-D)+/(R*D))**2 FOR RADIUS R(K) AND DIST-
!   ANCE D(K).  THUS
!
!        QX = (SWQX*SW - SWQ*SWX)/SW**2  AND
!        QY = (SWQY*SW - SWQ*SWY)/SW**2
!
!   WHERE SWQX AND SWX ARE PARTIAL DERIVATIVES WITH RESPECT
!   TO X OF SWQ AND SW, RESPECTIVELY.  SWQY AND SWY ARE DE-
!   FINED SIMILARLY.
!
            sw = 0.
            swx = 0.
            swy = 0.
            swq = 0.
            swqx = 0.
            swqy = 0.
!
! OUTER LOOP ON CELLS (I,J).
!
            do j = jmin , jmax
               do i = imin , imax
                  k = lcell(i,j)
                  if ( k==0 ) cycle
!
! INNER LOOP ON NODES K.
!
 5                delx = xp - x(k)
                  dely = yp - y(k)
                  dxsq = delx*delx
                  dysq = dely*dely
                  ds = dxsq + dysq
                  rs = rsq(k)
                  if ( ds<rs ) then
                     if ( ds==0. ) goto 20
                     rds = rs*ds
                     rd = sqrt(rds)
                     w = (rs+ds-rd-rd)/rds
                     t = 2.*(rd-rs)/(ds*rds)
                     wx = delx*t
                     wy = dely*t
                     qkx = 2.*a(1,k)*delx + a(2,k)*dely
                     qky = a(2,k)*delx + 2.*a(3,k)*dely
                     qk = (qkx*delx+qky*dely)/2.
                     qkx = qkx + a(4,k)
                     qky = qky + a(5,k)
                     qk = qk + a(4,k)*delx + a(5,k)*dely + f(k)
                     sw = sw + w
                     swx = swx + wx
                     swy = swy + wy
                     swq = swq + w*qk
                     swqx = swqx + wx*qk + w*qkx
                     swqy = swqy + wy*qk + w*qky
                  endif
!
! BOTTOM OF LOOP ON NODES IN CELL (I,J).
!
                  kp = k
                  k = lnext(kp)
                  if ( k/=kp ) goto 5
               enddo
            enddo
!
! SW = 0 IFF P IS NOT WITHIN THE RADIUS R(K) FOR ANY NODE K.
!
            if ( sw==0. ) goto 50
            q = swq/sw
            sws = sw*sw
            qx = (swqx*sw-swq*swx)/sws
            qy = (swqy*sw-swq*swy)/sws
            ier = 0
            return
!
! (PX,PY) = (X(K),Y(K))
!
 20         q = f(k)
            qx = a(4,k)
            qy = a(5,k)
            ier = 0
            return
         endif
!
! NO CELLS CONTAIN A POINT WITHIN RMAX OF P, OR
!   SW = 0 AND THUS DS .GE. RSQ(K) FOR ALL K.
!
 50      q = 0.
         qx = 0.
         qy = 0.
         ier = 2
      endif
      end subroutine qs2grd
!*==GETNP2.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine getnp2(px,py,x,y,nr,lcell,lnext,xmin,ymin,dx,dy,np,dsq)
      implicit none
!*--********************************************************************
!A INPUT  - PX
!A INPUT  - PY
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - NR
!A INPUT  - LCELL
!A OUTPUT - LNEXT
!A INPUT  - XMIN
!A INPUT  - YMIN
!A INPUT  - DX
!A INPUT  - DY
!A OUTPUT - NP
!A OUTPUT - DSQ
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   QSHEP2
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DELX     DELY     FIRST    I        I0       I1       I2       IMAX     IMIN     J        J0       J1       J2       JMAX     JMIN     L        LMIN     LN       R        RSMIN    RSQ      XP       YP
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real delx , dely , r , rsmin , rsq , xp , yp
      integer i , i0 , i1 , i2 , imax , imin , j , j0 , j1 , j2 , jmax ,&
     &        jmin , l , lmin , ln
!*** End of declarations inserted by SPAG
      integer nr , lcell(nr,nr) , lnext(*) , np
      real px , py , x(*) , y(*) , xmin , ymin , dx , dy , dsq
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!
!   GIVEN A SET OF N NODES AND THE DATA STRUCTURE DEFINED IN
! SUBROUTINE STORE2, THIS SUBROUTINE USES THE CELL METHOD TO
! FIND THE CLOSEST UNMARKED NODE NP TO A SPECIFIED POINT P.
! NP IS THEN MARKED BY SETTING LNEXT(NP) TO -LNEXT(NP).  (A
! NODE IS MARKED IF AND ONLY IF THE CORRESPONDING LNEXT ELE-
! MENT IS NEGATIVE.  THE ABSOLUTE VALUES OF LNEXT ELEMENTS,
! HOWEVER, MUST BE PRESERVED.)  THUS, THE CLOSEST M NODES TO
! P MAY BE DETERMINED BY A SEQUENCE OF M CALLS TO THIS ROU-
! TINE.  NOTE THAT IF THE NEAREST NEIGHBOR TO NODE K IS TO
! BE DETERMINED (PX = X(K) AND PY = Y(K)), THEN K SHOULD BE
! MARKED BEFORE THE CALL TO THIS ROUTINE.
!   THE SEARCH IS BEGUN IN THE CELL CONTAINING (OR CLOSEST
! TO) P AND PROCEEDS OUTWARD IN RECTANGULAR LAYERS UNTIL ALL
! CELLS WHICH CONTAIN POINTS WITHIN DISTANCE R OF P HAVE
! BEEN SEARCHED, WHERE R IS THE DISTANCE FROM P TO THE FIRST
! UNMARKED NODE ENCOUNTERED (INFINITE IF NO UNMARKED NODES
! ARE PRESENT).
!
! ON INPUT --
!
!       PX,PY = CARTESIAN COORDINATES OF THE POINT P WHOSE
!               NEAREST UNMARKED NEIGHBOR IS TO BE FOUND.
!
!       X,Y = ARRAYS OF LENGTH N, FOR N .GE. 2, CONTAINING
!             THE CARTESIAN COORDINATES OF THE NODES.
!
!       NR = NUMBER OF ROWS AND COLUMNS IN THE CELL GRID.
!            NR .GE. 1.
!
!       LCELL = NR BY NR ARRAY OF NODAL INDICES ASSOCIATED
!               WITH CELLS.
!
!       LNEXT = ARRAY OF LENGTH N CONTAINING NEXT-NODE INDI-
!               CES (OR THEIR NEGATIVES).
!
!       XMIN,YMIN,DX,DY = MINIMUM NODAL COORDINATES AND CELL
!                         DIMENSIONS.  DX AND DY MUST BE
!                         POSITIVE.
!
!   INPUT PARAMETERS OTHER THAN LNEXT ARE NOT ALTERED BY
! THIS ROUTINE.  WITH THE EXCEPTION OF (PX,PY) AND THE SIGNS
! OF LNEXT ELEMENTS, THESE PARAMETERS SHOULD BE UNALTERED
! FROM THEIR VALUES ON OUTPUT FROM SUBROUTINE STORE2.
!
! ON OUTPUT --
!
!       NP = INDEX (FOR X AND Y) OF THE NEAREST UNMARKED
!            NODE TO P, OR 0 IF ALL NODES ARE MARKED OR NR
!            .LT. 1 OR DX .LE. 0 OR DY .LE. 0.  LNEXT(NP)
!            .LT. 0 IF NP .NE. 0.
!
!       DSQ = SQUARED EUCLIDEAN DISTANCE BETWEEN P AND NODE
!             NP, OR 0 IF NP = 0.
!
! MODULES REQUIRED BY GETNP2 -- NONE
!
! INTRINSIC FUNCTIONS CALLED BY GETNP2 -- IABS, IFIX, SQRT
!
!***********************************************************
!
      logical first
      xp = px
      yp = py
!
! TEST FOR INVALID INPUT PARAMETERS.
!
      if ( nr<1 .or. dx<=0. .or. dy<=0. ) goto 200
!
! INITIALIZE PARAMETERS --
!
!   FIRST = TRUE IFF THE FIRST UNMARKED NODE HAS YET TO BE
!           ENCOUNTERED,
!   IMIN,IMAX,JMIN,JMAX = CELL INDICES DEFINING THE RANGE OF
!                         THE SEARCH,
!   DELX,DELY = PX-XMIN AND PY-YMIN,
!   I0,J0 = CELL CONTAINING OR CLOSEST TO P,
!   I1,I2,J1,J2 = CELL INDICES OF THE LAYER WHOSE INTERSEC-
!                 TION WITH THE RANGE DEFINED BY IMIN,...,
!                 JMAX IS CURRENTLY BEING SEARCHED.
!
      first = .true.
      imin = 1
      imax = nr
      jmin = 1
      jmax = nr
      delx = xp - xmin
      dely = yp - ymin
      i0 = ifix(delx/dx) + 1
      if ( i0<1 ) i0 = 1
      if ( i0>nr ) i0 = nr
      j0 = ifix(dely/dy) + 1
      if ( j0<1 ) j0 = 1
      if ( j0>nr ) j0 = nr
      i1 = i0
      i2 = i0
      j1 = j0
      j2 = j0
!
! OUTER LOOP ON LAYERS, INNER LOOP ON LAYER CELLS, EXCLUDING
!   THOSE OUTSIDE THE RANGE (IMIN,IMAX) X (JMIN,JMAX).
!
 100  do j = j1 , j2
         if ( j>jmax ) exit
         if ( j>=jmin ) then
            do i = i1 , i2
               if ( i>imax ) exit
               if ( i<imin ) cycle
               if ( j/=j1 .and. j/=j2 .and. i/=i1 .and. i/=i2 ) cycle
!
! SEARCH CELL (I,J) FOR UNMARKED NODES L.
!
               l = lcell(i,j)
               if ( l==0 ) cycle
!
!   LOOP ON NODES IN CELL (I,J).
!
 110           ln = lnext(l)
               if ( ln>=0 ) then
!
!   NODE L IS NOT MARKED.
!
                  rsq = (x(l)-xp)**2 + (y(l)-yp)**2
                  if ( first ) then
!
!   NODE L IS THE FIRST UNMARKED NEIGHBOR OF P ENCOUNTERED.
!     INITIALIZE LMIN TO THE CURRENT CANDIDATE FOR NP, AND
!     RSMIN TO THE SQUARED DISTANCE FROM P TO LMIN.  IMIN,
!     IMAX, JMIN, AND JMAX ARE UPDATED TO DEFINE THE SMAL-
!     LEST RECTANGLE CONTAINING A CIRCLE OF RADIUS R =
!     SQRT(RSMIN) CENTERED AT P, AND CONTAINED IN (1,NR) X
!     (1,NR) (EXCEPT THAT, IF P IS OUTSIDE THE RECTANGLE
!     DEFINED BY THE NODES, IT IS POSSIBLE THAT IMIN .GT.
!     NR, IMAX .LT. 1, JMIN .GT. NR, OR JMAX .LT. 1).  FIRST
!     IS RESET TO FALSE.
!
                     lmin = l
                     rsmin = rsq
                     r = sqrt(rsmin)
                     imin = ifix((delx-r)/dx) + 1
                     if ( imin<1 ) imin = 1
                     imax = ifix((delx+r)/dx) + 1
                     if ( imax>nr ) imax = nr
                     jmin = ifix((dely-r)/dy) + 1
                     if ( jmin<1 ) jmin = 1
                     jmax = ifix((dely+r)/dy) + 1
                     if ( jmax>nr ) jmax = nr
                     first = .false.
!
!   TEST FOR NODE L CLOSER THAN LMIN TO P.
!
                  elseif ( rsq<rsmin ) then
!
!   UPDATE LMIN AND RSMIN.
!
                     lmin = l
                     rsmin = rsq
                  endif
               endif
!
!   TEST FOR TERMINATION OF LOOP ON NODES IN CELL (I,J).
!
               if ( iabs(ln)/=l ) then
                  l = iabs(ln)
                  goto 110
               endif
            enddo
         endif
      enddo
!
! TEST FOR TERMINATION OF LOOP ON CELL LAYERS.
!
      if ( i1>imin .or. i2<imax .or. j1>jmin .or. j2<jmax ) then
         i1 = i1 - 1
         i2 = i2 + 1
         j1 = j1 - 1
         j2 = j2 + 1
         goto 100
!
! UNLESS NO UNMARKED NODES WERE ENCOUNTERED, LMIN IS THE
!   CLOSEST UNMARKED NODE TO P.
!
      elseif ( .not.(first) ) then
         np = lmin
         dsq = rsmin
         lnext(lmin) = -lnext(lmin)
         return
      endif
!
! ERROR -- NR, DX, OR DY IS INVALID OR ALL NODES ARE MARKED.
!
 200  np = 0
      dsq = 0.
      end subroutine getnp2
!*==GIVENS.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine givens(a,b,c,s)
      implicit none
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
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
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
! ON INPUT --
!
!       A,B = COMPONENTS OF THE 2-VECTOR TO BE ROTATED.
!
! ON OUTPUT --
!
!       A = VALUE OVERWRITTEN BY R = +/-SQRT(A*A + B*B)
!
!       B = VALUE OVERWRITTEN BY A VALUE Z WHICH ALLOWS C
!           AND S TO BE RECOVERED AS FOLLOWS --
!             C = SQRT(1-Z*Z), S=Z     IF ABS(Z) .LE. 1.
!             C = 1/Z, S = SQRT(1-C*C) IF ABS(Z) .GT. 1.
!
!       C = +/-(A/R)
!
!       S = +/-(B/R)
!
! MODULES REQUIRED BY GIVENS -- NONE
!
! INTRINSIC FUNCTIONS CALLED BY GIVENS - ABS, SQRT
!
!***********************************************************
!
      real aa , bb , r , u , v
!
! LOCAL PARAMETERS --
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
!   SIGN(A)*SIGN(B).
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
! STORE R IN A.
!
      a = sqrt(.25+v*v)*u
      s = bb/a
      c = v*(s+s)
!
! NOTE THAT R HAS THE SIGN OF B, S .GT. 0, AND C HAS
!   SIGN(A)*SIGN(B).
!
      b = 1.
      if ( c/=0. ) b = 1./c
      return
99999 end subroutine givens
!*==ROTATE.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine rotate(n,c,s,x,y)
      implicit none
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
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!
!                                            ( C  S)
!   THIS ROUTINE APPLIES THE GIVENS ROTATION (     ) TO THE
!                                            (-S  C)
!               (X(1) ... X(N))
! 2 BY N MATRIX (             ).
!               (Y(1) ... Y(N))
!
! ON INPUT --
!
!       N = NUMBER OF COLUMNS TO BE ROTATED.
!
!       C,S = ELEMENTS OF THE GIVENS ROTATION.  THESE MAY BE
!             DETERMINED BY SUBROUTINE GIVENS.
!
!       X,Y = ARRAYS OF LENGTH .GE. N CONTAINING THE VECTORS
!             TO BE ROTATED.
!
! PARAMETERS N, C, AND S ARE NOT ALTERED BY THIS ROUTINE.
!
! ON OUTPUT --
!
!       X,Y = ROTATED VECTORS.
!
! MODULES REQUIRED BY ROTATE -- NONE
!
!***********************************************************
!
      integer i
      real xi , yi
!
! LOCAL PARAMETERS --
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
!*==SETUP2.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine setup2(xk,yk,fk,xi,yi,fi,s1,s2,r,row)
      implicit none
!*--********************************************************************
!A INPUT  - XK
!A INPUT  - YK
!A INPUT  - FK
!A INPUT  - XI
!A INPUT  - YI
!A INPUT  - FI
!A INPUT  - S1
!A INPUT  - S2
!A INPUT  - R
!A OUTPUT - ROW
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   QSHEP2
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  D        DX       DXSQ     DY       DYSQ     I        W        W1       W2
! uses PARAMs *** NONE ****
!*++********************************************************************
      real xk , yk , fk , xi , yi , fi , s1 , s2 , r , row(6)
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!
!   THIS ROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED RE-
! GRESSION MATRIX FOR A WEIGHTED LEAST-SQUARES FIT OF A
! QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES F, WHERE
! Q(XK,YK) = FK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS) ARE
! SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LINEAR
! TERMS) ARE SCALED BY 1/S1.  THE WEIGHT IS (R-D)/(R*D) IF
! R .GT. D AND 0 IF R .LE. D, WHERE D IS THE DISTANCE
! BETWEEN NODES I AND K.
!
! ON INPUT --
!
!       XK,YK,FK = COORDINATES AND DATA VALUE AT NODE K --
!                  INTERPOLATED BY Q.
!
!       XI,YI,FI = COORDINATES AND DATA VALUE AT NODE I.
!
!       S1,S2 = RECIPROCALS OF THE SCALE FACTORS.
!
!       R = RADIUS OF INFLUENCE ABOUT NODE K DEFINING THE
!           WEIGHT.
!
!       ROW = ARRAY OF LENGTH 6.
!
! INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
! ON OUTPUT --
!
!       ROW = VECTOR CONTAINING A ROW OF THE AUGMENTED
!             REGRESSION MATRIX.
!
! MODULES REQUIRED BY SETUP2 -- NONE
!
! INTRINSIC FUNCTION CALLED BY SETUP2 -- SQRT
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
      row(6) = (fi-fk)*w
      return
99999 end subroutine setup2
!*==STORE2.f90  processed by SPAG 7.51DB at 21:35 on 14 Mar 2022
      subroutine store2(n,x,y,nr,lcell,lnext,xmin,ymin,dx,dy,ier)
      implicit none
!*--********************************************************************
!A INPUT  - N
!A INPUT  - X
!A INPUT  - Y
!A INPUT  - NR
!A OUTPUT - LCELL
!A OUTPUT - LNEXT
!A OUTPUT - XMIN
!A OUTPUT - YMIN
!A OUTPUT - DX
!A OUTPUT - DY
!A OUTPUT - IER
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calls       ** NOTHING **
! called by   QSHEP2
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  DELX     DELY     I        J        K        KB       L        NN       NNR      NP1      XMN      XMX      YMN      YMX
! uses PARAMs *** NONE ****
!*++********************************************************************
!*** Start of declarations inserted by SPAG
      real delx , dely , xmn , xmx , ymn , ymx
      integer i , j , k , kb , l , nn , nnr , np1
!*** End of declarations inserted by SPAG
      integer n , nr , lcell(nr,nr) , lnext(n) , ier
      real x(n) , y(n) , xmin , ymin , dx , dy
!
!***********************************************************
!
!                                               ROBERT RENKA
!                                       UNIV. OF NORTH TEXAS
!                                             (817) 565-2767
!
!   GIVEN A SET OF N ARBITRARILY DISTRIBUTED NODES IN THE
! PLANE, THIS SUBROUTINE CREATES A DATA STRUCTURE FOR A
! CELL-BASED METHOD OF SOLVING CLOSEST-POINT PROBLEMS.  THE
! SMALLEST RECTANGLE CONTAINING THE NODES IS PARTITIONED
! INTO AN NR BY NR UNIFORM GRID OF CELLS, AND NODES ARE AS-
! SOCIATED WITH CELLS.  IN PARTICULAR, THE DATA STRUCTURE
! STORES THE INDICES OF THE NODES CONTAINED IN EACH CELL.
! FOR A UNIFORM RANDOM DISTRIBUTION OF NODES, THE NEAREST
! NODE TO AN ARBITRARY POINT CAN BE DETERMINED IN CONSTANT
! EXPECTED TIME.
!
! ON INPUT --
!
!       N = NUMBER OF NODES.  N .GE. 2.
!
!       X,Y = ARRAYS OF LENGTH N CONTAINING THE CARTESIAN
!             COORDINATES OF THE NODES.
!
!       NR = NUMBER OF ROWS AND COLUMNS IN THE GRID.  THE
!            CELL DENSITY (AVERAGE NUMBER OF NODES PER CELL)
!            IS D = N/(NR**2).  A RECOMMENDED VALUE, BASED
!            ON EMPIRICAL EVIDENCE, IS D = 3 -- NR =
!            SQRT(N/3).  NR .GE. 1.
!
! THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.
!
!       LCELL = ARRAY OF LENGTH .GE. NR**2.
!
!       LNEXT = ARRAY OF LENGTH .GE. N.
!
! ON OUTPUT --
!
!       LCELL = NR BY NR CELL ARRAY SUCH THAT LCELL(I,J)
!               CONTAINS THE INDEX (FOR X AND Y) OF THE
!               FIRST NODE (NODE WITH SMALLEST INDEX) IN
!               CELL (I,J), OR LCELL(I,J) = 0 IF NO NODES
!               ARE CONTAINED IN THE CELL.  THE UPPER RIGHT
!               CORNER OF CELL (I,J) HAS COORDINATES (XMIN+
!               I*DX,YMIN+J*DY).  LCELL IS NOT DEFINED IF
!               IER .NE. 0.
!
!       LNEXT = ARRAY OF NEXT-NODE INDICES SUCH THAT
!               LNEXT(K) CONTAINS THE INDEX OF THE NEXT NODE
!               IN THE CELL WHICH CONTAINS NODE K, OR
!               LNEXT(K) = K IF K IS THE LAST NODE IN THE
!               CELL FOR K = 1,...,N.  (THE NODES CONTAINED
!               IN A CELL ARE ORDERED BY THEIR INDICES.)
!               IF, FOR EXAMPLE, CELL (I,J) CONTAINS NODES
!               2, 3, AND 5 (AND NO OTHERS), THEN LCELL(I,J)
!               = 2, LNEXT(2) = 3, LNEXT(3) = 5, AND
!               LNEXT(5) = 5.  LNEXT IS NOT DEFINED IF
!               IER .NE. 0.
!
!       XMIN,YMIN = CARTESIAN COORDINATES OF THE LOWER LEFT
!                   CORNER OF THE RECTANGLE DEFINED BY THE
!                   NODES (SMALLEST NODAL COORDINATES) UN-
!                   LESS IER = 1.  THE UPPER RIGHT CORNER IS
!                   (XMAX,YMAX) FOR XMAX = XMIN + NR*DX AND
!                   YMAX = YMIN + NR*DY.
!
!       DX,DY = DIMENSIONS OF THE CELLS UNLESS IER = 1.  DX
!               = (XMAX-XMIN)/NR AND DY = (YMAX-YMIN)/NR
!               WHERE XMIN, XMAX, YMIN, AND YMAX ARE THE
!               EXTREMA OF X AND Y.
!
!       IER = ERROR INDICATOR --
!             IER = 0 IF NO ERRORS WERE ENCOUNTERED.
!             IER = 1 IF N .LT. 2 OR NR .LT. 1.
!             IER = 2 IF DX = 0 OR DY = 0.
!
! MODULES REQUIRED BY STORE2 -- NONE
!
! INTRINSIC FUNCTIONS CALLED BY STORE2 -- FLOAT, IFIX
!
!***********************************************************
!
      nn = n
      nnr = nr
      if ( nn<2 .or. nnr<1 ) then
!
! INVALID INPUT PARAMETER
!
         ier = 1
         return
      else
!
! COMPUTE THE DIMENSIONS OF THE RECTANGLE CONTAINING THE
!   NODES.
!
         xmn = x(1)
         xmx = xmn
         ymn = y(1)
         ymx = ymn
         do k = 2 , nn
            if ( x(k)<xmn ) xmn = x(k)
            if ( x(k)>xmx ) xmx = x(k)
            if ( y(k)<ymn ) ymn = y(k)
            if ( y(k)>ymx ) ymx = y(k)
         enddo
         xmin = xmn
         ymin = ymn
!
! COMPUTE CELL DIMENSIONS AND TEST FOR ZERO AREA.
!
         delx = (xmx-xmn)/float(nnr)
         dely = (ymx-ymn)/float(nnr)
         dx = delx
         dy = dely
         if ( delx==0. .or. dely==0. ) then
!
! DX = 0 OR DY = 0
!
            ier = 2
            goto 99999
         endif
      endif
!
! INITIALIZE LCELL.
!
      do j = 1 , nnr
         do i = 1 , nnr
            lcell(i,j) = 0
         enddo
      enddo
!
! LOOP ON NODES, STORING INDICES IN LCELL AND LNEXT.
!
      np1 = nn + 1
      do k = 1 , nn
         kb = np1 - k
         i = ifix((x(kb)-xmn)/delx) + 1
         if ( i>nnr ) i = nnr
         j = ifix((y(kb)-ymn)/dely) + 1
         if ( j>nnr ) j = nnr
         l = lcell(i,j)
         lnext(kb) = l
         if ( l==0 ) lnext(kb) = kb
         lcell(i,j) = kb
      enddo
!
! NO ERRORS ENCOUNTERED
!
      ier = 0
      return
99999 end subroutine store2
end module alg660
