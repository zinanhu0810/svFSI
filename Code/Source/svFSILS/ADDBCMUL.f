!--------------------------------------------------------------------
!     Created by Mahdi Esmaily Moghadam
!     contact memt63@gmail.com for reporting the bugs.
!--------------------------------------------------------------------
!
!     UC Copyright Notice
!     This software is Copyright ©2012 The Regents of the University of
!     California. All Rights Reserved.
!
!     Permission to copy and modify this software and its documentation
!     for educational, research and non-profit purposes, without fee,
!     and without a written agreement is hereby granted, provided that
!     the above copyright notice, this paragraph and the following three
!     paragraphs appear in all copies.
!
!     Permission to make commercial use of this software may be obtained
!     by contacting:
!     Technology Transfer Office
!     9500 Gilman Drive, Mail Code 0910
!     University of California
!     La Jolla, CA 92093-0910
!     (858) 534-5815
!     invent@ucsd.edu
!
!     This software program and documentation are copyrighted by The
!     Regents of the University of California. The software program and
!     documentation are supplied "as is", without any accompanying
!     services from The Regents. The Regents does not warrant that the
!     operation of the program will be uninterrupted or error-free. The
!     end-user understands that the program was developed for research
!     purposes and is advised not to rely exclusively on the program for
!     any reason.
!
!     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
!     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
!     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
!     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
!     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
!     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
!     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
!     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
!     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
!
!--------------------------------------------------------------------
!     The contribution of coupled BCs is added to the matrix-vector
!     product operation. Depending on the type of operation (adding the
!     contribution or computing the PC contribution) different
!     coefficients are used.
!
!     AB: The resistance is stored in lhs%face(faIn)%res.
!     The current matrix vector product is in Y. The vector to be
!     multiplied and added is in X. The matrix is represented by res.
!     See for reference
!        - Moghadam et al. 2013 eq. 27
!          (https://doi.org/10.1016/j.jcp.2012.07.035)
!        - Moghadam et al. 2013b
!          (https://doi.org/10.1007/s00466-013-0868-1)
!
!--------------------------------------------------------------------

      SUBROUTINE ADDBCMUL(lhs, op_Type, dof, X, Y)
      INCLUDE "FSILS_STD.h"
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: op_type, dof
      REAL(KIND=LSRP), INTENT(IN) :: X(dof, lhs%nNo)
      REAL(KIND=LSRP), INTENT(INOUT) :: Y(dof, lhs%nNo)

      INTEGER(KIND=LSIP) faIn, faInCap, i, a, Ac, nsd
      REAL(KIND=LSRP) S, FSILS_DOTV

      REAL(KIND=LSRP), ALLOCATABLE :: v(:,:), vcap(:,:), coef(:)

      ALLOCATE(coef(lhs%nFaces), v(dof,lhs%nNo), vcap(dof,lhs%nNo))

!     Setting coef depending on adding resistance to stiffness or
!     computing preconditioner
      IF (op_Type .EQ. BCOP_TYPE_ADD) THEN
         coef = lhs%face%res
      ELSE IF(op_Type .EQ. BCOP_TYPE_PRE) THEN
         coef = -lhs%face%res/(1._LSRP + (lhs%face%res*lhs%face%nS))
      ELSE
         PRINT *, "FSILS: op_Type is not defined"
         STOP "FSILS: FATAL ERROR"
      END IF

      DO faIn=1, lhs%nFaces

!        If the face is virtual, don't add anything to tangent matrix.
         IF (lhs%face(faIn)%vrtual) CYCLE

!        In the following calculations, we are computing the product of
!        the coupled BC tangent contribution with the vector X (refer to
!        Moghadam et al. 2013 eq. 27). This is computed by first
!        constructing the vector v, which one of the integrals found in
!        the expression, int{N_A * n_i} dGamma. Then, v is dotted with X
!        to yield a quantity S. Then S is multiplied by v again, and
!        also multiplied by the appropriate coefficients in the
!        expression. The calculations are complicated somewhat if there
!        is a capping surface, but these complications are explained
!        below.

!        Calculating S, which is the inner product of the right integral
!        (v) and the vector to be multiplied (X).
         nsd = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%coupledFlag) THEN
            IF (lhs%face(faIn)%sharedFlag) THEN
               v = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     v(i,Ac) = lhs%face(faIn)%valM(i,a)
                  END DO
               END DO
               S = coef(faIn)*FSILS_DOTV(dof,lhs%mynNo, lhs%commu, v, X)

!              If a virtual face caps this face, add contribution to S
!              Explanation: If the coupled surface is virtually capped
!              to compute flow rate then the right integral should be
!              over the capped surface + capping surface, while the left
!              integral should be over the uncapped surface (because we
!              do not apply a pressure to the capping surface). We can
!              achieve this by adding the capping face's contribution
               IF (lhs%face(faIn)%faInCap .NE. 0) THEN
                     faInCap = lhs%face(faIn)%faInCap
                     IF(.NOT. lhs%face(faInCap)%coupledFlag) THEN
                        PRINT*, 'ADDBCMUL(): Cap face is not coupled.',
     2                     'Probably cap face has zero resistance.'
                        STOP "FSILS: FATAL ERROR"
                     END IF
                     vcap = 0._LSRP
                     DO a=1, lhs%face(faInCap)%nNo
                        Ac = lhs%face(faInCap)%glob(a)
                        DO i=1, nsd
                           vcap(i,Ac) = lhs%face(faInCap)%valM(i,a)
                        END DO
                     END DO
                     S = S + coef(faIn)*FSILS_DOTV(dof,lhs%mynNo,
     2                  lhs%commu, vcap, X)
               END IF

!              Add S times left integral to the current matrix-vector
!              product Y. Note that we do not add the capping surface's
!              contribution to v (vcap), since the left integral should
!              be over the uncapped surface.
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S
                  END DO
               END DO
            ELSE
               S = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     S = S + lhs%face(faIn)%valM(i,a)*X(i,Ac)
                  END DO
               END DO
               IF (lhs%face(faIn)%faInCap .NE. 0) THEN
                     faInCap = lhs%face(faIn)%faInCap
                     IF(.NOT. lhs%face(faInCap)%coupledFlag) THEN
                        PRINT*, 'ADDBCMUL(): Cap face is not coupled.',
     2                     'Probably cap face has zero resistance.'
                     STOP "FSILS: FATAL ERROR"
                     END IF

                     DO a=1, lhs%face(faInCap)%nNo
                        Ac = lhs%face(faInCap)%glob(a)
                        DO i=1, nsd
                           S = S + lhs%face(faInCap)%valM(i,a)*X(i,Ac)
                        END DO
                     END DO
               END IF

!              Multiply S by the resistance or related quantity if
!              preconditioning
               S = coef(faIn)*S

!              Add S times left integral to the current matrix-vector
!              product Y. Note that we do not add the capping surface's
!              contribution to v (vcap), since the left integral should
!              be over the uncapped surface.
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     Y(i,Ac) = Y(i,Ac) + lhs%face(faIn)%valM(i,a)*S
                  END DO
               END DO
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE ADDBCMUL
!####################################################################
