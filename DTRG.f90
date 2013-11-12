!INCLUDE 'MATHIO.F90'
!INCLUDE 'TENSOR.F90'
! ############### CONST ###################
MODULE CONST
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL
	USE CONST
	REAL    :: THETA = 1.*PI ! theta-term
	REAL    :: CROSS = 0.    ! crossing: 1. = allow, 0. = avoid
	REAL    :: BETA = 0.    ! inverse temperature
	INTEGER :: LEN = 4 
	INTEGER :: MAX_CUT = 6  ! 16
	REAL    :: MAX_ERR = 0.
	INTEGER :: SWEEPS = 1
END MODULE MODEL
! ################ MATH ####################
MODULE MATH
	USE TENSORIAL
	TYPE HUGE_TENSOR
		REAL :: LEV
		TYPE(TENSOR) :: TEN
	END TYPE HUGE_TENSOR
CONTAINS
! to huge tensor
FUNCTION TOHUGE(T) RESULT (M)
	TYPE(TENSOR), INTENT(IN) :: T
	TYPE(HUGE_TENSOR) :: M
	
	M%LEV = 0.
	M%TEN = T
	CALL NORMALIZE(M)
END FUNCTION TOHUGE
! normalize a huge tensor
SUBROUTINE NORMALIZE(M)
	TYPE(HUGE_TENSOR), INTENT(INOUT) :: M
	! local variables
	REAL :: R
	
	R = SUM(ABS(M%TEN%VALS)) ! get ratio
	IF (R >= 1.E-12) THEN
		M%TEN%VALS = M%TEN%VALS/R ! rescale
		M%LEV = M%LEV + LOG(R) ! update level
	END IF
END SUBROUTINE NORMALIZE
! end of module MATH
END MODULE MATH
! ############## PHYSICS ###################
MODULE PHYSICS
	USE MATH
! tensor conventions
! MPO (by default B-type):
!   A-type       B-type
!     ╭            ╮
!     3            3
! ─ 2 A 1 ─    ─ 1 B 2 ─
!     4            4
!     ╰            ╯
! MPS (by default A-type):
!   A-type       B-type
!     │            │
!     2            2
! ─ 1 A 3 ─    ─ 3 B 1 ─
CONTAINS
! ------------ set MPO tensor ---------------
! square lattice MPO
SUBROUTINE SET_MPO(T)
! output: T - MPO tensor (by default B-type)
!     │
!     3
! ─ 1 T 2 ─
!     4
!     │ 
! leg-1,2: physical, leg-3,4: internal
	USE MODEL
	TYPE(TENSOR), INTENT(OUT) :: T ! MPO tensor output
	! local variables
	TYPE(TENSOR) :: X, Y, UA, UB, S
	COMPLEX, ALLOCATABLE :: WA(:,:)
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA * ZI
	B = BETA * Z1
	X =  TENSOR([2,2,2,2],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[EXP(2*B),EXP(Q/4.),EXP(Q/4.),EXP(Z0),EXP(Q/4.),CROSS*EXP(-2*B - Q/2.),EXP(Z0),EXP(-Q/4.),EXP(Q/4.),EXP(Z0),CROSS*EXP(-2*B + Q/2.),EXP(-Q/4.),EXP(Z0),EXP(-Q/4.),EXP(-Q/4.),EXP(2*B)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! SVD of X to UA - S - UB
	CALL SVD(X,[1,4],[3,2],UA,UB,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	! attach the half S to U
	UA = TEN_PROD(S,UA,[2],[3])
	UB = TEN_PROD(S,UB,[2],[3])
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	UA = TEN_PROD(UA,Y,[3],[3])
	UB = TEN_PROD(UB,Y,[3],[3])
	T = TEN_TRANS(TEN_PROD(UB,UA,[2,4],[4,2]),[2,4,1,3])
END SUBROUTINE SET_MPO
! ----------- DMRG -----------
! DMRG controlling routine
SUBROUTINE DMRG(TB, M)
! input: T - MPO site tensor
! output: WA - MPS tensor
	USE MODEL
	TYPE(TENSOR), INTENT(IN)  :: TB
	TYPE(HUGE_TENSOR), INTENT(OUT) :: M
	! local tensor
	TYPE(TENSOR) :: TA, WA(LEN), WB(LEN), DA(LEN), DB(LEN), P
	! local variables
	INTEGER :: L, ITER
	COMPLEX :: F, FA(LEN), FB(LEN)
	
	! check validity of system size LEN
	IF (MODULO(LEN,2) == 1 .OR. LEN < 4) THEN
		WRITE (*,'(A)') 'DMRG::xlen: LEN must be even and greater than 4.'
		STOP
	END IF
	! set TA from TB MPO
	TA = TEN_TRANS(TB,[2,1,3,4]) ! given by 1 <-> 2
	! iDMRG (warm up)
	DO L = 1, LEN/2
		CALL IDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
		CALL SHOW_DMRG_STATUS('I',L,F,P)
	END DO
	! fDMRG (first sweep)
	DO L = LEN/2+1, LEN
		CALL FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
		CALL SHOW_DMRG_STATUS('F',L,F,P)
	END DO
	! fDMRG sweeps
	DO ITER = 1, SWEEPS
		! fFDMRG (backward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, TB, TA, DB, DA, FB, FA, WB, WA, P, F)
			CALL SHOW_DMRG_STATUS('B',L,F,P)
		END DO
		! fFDMRG (forward sweep)
		DO L = 1, LEN
			CALL FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
			CALL SHOW_DMRG_STATUS('F',L,F,P)
		END DO
	END DO
	! fDMRG (last sweep)
	DO L = 1, LEN/2 - 1
		CALL FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
		CALL SHOW_DMRG_STATUS('F',L,F,P)
	END DO
	CALL LDMRG(LEN/2, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
END SUBROUTINE DMRG
! infinite-size DMRG step
SUBROUTINE IDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
! perform one step of iDMRG at L
! input: L - lattice position
!        TA,TB - MPO
!        DA, DB - block MPO, FA, FB - level of block MPO
!        WA, WB - MPS, P - Schmidt spectrum 
!        F - free energy
! inout: DA, DB, WA, WB, S0, S will be updated
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: TA, TB
	TYPE(TENSOR), INTENT(INOUT) :: DA(:), DB(:), WA(:), WB(:), P
	COMPLEX, INTENT(INOUT) :: F, FA(:), FB(:)
	! local tensors
	TYPE(TENSOR) :: W, TS
	TYPE(TENSOR), SAVE :: P0
	! local variables
	INTEGER :: DPHY
	COMPLEX :: TVAL
	
	IF (L == 1) THEN
		! initialize tensors
		DPHY = TB%DIMS(3) ! get physical dim
		! set initial Schmidt spectrum
		P0 = EYE_TEN([1,1])
		P = EYE_TEN([DPHY,DPHY],SQRT(Z1/DPHY))
		! set initial MPS (boundary)
		!   2       2
		! 1 WA 3  3 WB 1
		! initial MPS transfers d.o.f. from leg 2 to 3, with leg 1 dummy
		WA(1) = TEN_PROD(TENSOR([1],[0],[Z1]),EYE_TEN([DPHY,DPHY]))
		WB(1) = WA(1) ! WB is the same as WA
		F = Z0 ! set an initial F
		! set initial block MPO (and rescale)
		DA(1) = TA
		DB(1) = TB
		CALL RESCALE(DA(1), Z0, FA(1))
		CALL RESCALE(DB(1), Z0, FB(1))
	ELSE ! L >= 2
		! estimate trial W
		P0%VALS = Z1/P0%VALS ! cal P0^(-1) -> P0
		W = MAKE_W5(WA(L-1), WB(L-1), P, P0)
		P0 = P ! P0 is used, update to P
		TS = MAKE_TS(TA, TB, DA(L-1), DB(L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(L-1)+FB(L-1)+LOG(TVAL) ! calculate F
		! SVD split W, and update WA, WB, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(L),P,MAX_CUT,MAX_ERR)
		! update DA, DB and rescale
		DA(L) = NEW_DX(TA, WA(L), DA(L-1))
		DB(L) = NEW_DX(TB, WB(L), DB(L-1))
		CALL RESCALE(DA(L), FA(L-1), FA(L))
		CALL RESCALE(DB(L), FB(L-1), FB(L))
	END IF
END SUBROUTINE IDMRG
! finite-size DMRG step (forward convention)
SUBROUTINE FDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
! perform one step of fDMRG at L (forward dirction)
! input: L - lattice position
!        TA,TB - MPO
!        DA, DB - block MPO, LA, LB - level of block MPO
!        WA, WB - MPS, P - Schmidt spectrum 
!        F - free energy
! inout: DA, WA, S will be updated
! on output: DB is not touched, WB is destroyed
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: TA, TB
	TYPE(TENSOR), INTENT(INOUT) :: DA(:), DB(:), WA(:), WB(:), P
	COMPLEX, INTENT(INOUT) :: F, FA(:), FB(:)
	! local tensors
	TYPE(TENSOR) :: W, TS
	COMPLEX :: TVAL
	COMPLEX, SAVE :: F0

	IF (L == 1) THEN
		! WA(1), no update, obtained from the last sweep
		! update DA, to restart
		! SSB treatment to be implemented here
		DA(1) = NEW_DX(TA, WA(1))
		CALL RESCALE(DA(1), Z0, FA(1)) ! rescale
		F = F0 ! retrieve F from last save
	ELSEIF (L == LEN-1) THEN
		! estimate trial W
		W = MAKE_W3(P, WB(2), WB(1))
		TS = MAKE_TS(TA, TB, DA(LEN-2)) ! construct TS (boundary)
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(LEN-2)+LOG(TVAL) ! calculate F
		F0 = F ! save F at the last site
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(LEN-1),WB(1),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(LEN-1) = NEW_DX(TA, WA(LEN-1), DA(LEN-2))
		CALL RESCALE(DA(LEN-1), FA(LEN-2), FA(LEN-1))
	ELSEIF (L == LEN) THEN
		! update the ending WA by P*WB
		WA(LEN) = TEN_TRANS(TEN_PROD(P,WB(1),[2],[3]),[1,3,2])
		! update DA and rescale
		DA(LEN) = NEW_DX(TA, WA(LEN), DA(LEN-1))
		CALL RESCALE(DA(LEN), FA(LEN-1), FA(LEN))
		F = FA(LEN) ! return F for the whole lattice
	ELSE ! 2 <= L <= LEN-2
		! estimate trial W
		W = MAKE_W3(P, WB(LEN-L+1), WB(LEN-L))
		TS = MAKE_TS(TA, TB, DA(L-1), DB(LEN-L-1)) ! construct TS
		TVAL = ANNEAL(TS, W) ! anneal W by TS
		F = FA(L-1)+FB(LEN-L-1)+LOG(TVAL) ! calculate F
		! SVD split W, update WA, P
		CALL SVD(W,[1,2],[3,4],WA(L),WB(LEN-L),P,MAX_CUT,MAX_ERR)
		! update DA and rescale
		DA(L) = NEW_DX(TA, WA(L), DA(L-1))
		CALL RESCALE(DA(L), FA(L-1), FA(L))
	END IF
END SUBROUTINE FDMRG
! last DMRG step
SUBROUTINE LDMRG(L, TA, TB, DA, DB, FA, FB, WA, WB, P, F)
! perform one step of fDMRG at L (forward dirction)
! input: L - lattice position
!        TA,TB - MPO
!        DA, DB - block MPO, LA, LB - level of block MPO
!        WA, WB - MPS, P - Schmidt spectrum 
!        F - free energy
	USE MODEL
	INTEGER, INTENT(IN)      :: L
	TYPE(TENSOR), INTENT(IN) :: TA, TB
	TYPE(TENSOR), INTENT(INOUT) :: DA(:), DB(:), WA(:), WB(:), P
	COMPLEX, INTENT(INOUT) :: F, FA(:), FB(:)
	! local tensors
	TYPE(TENSOR) :: W, TS, U
	COMPLEX :: TVAL
	
	! estimate trial W
	W = MAKE_W3(P, WB(LEN-L+1), WB(LEN-L))
	TS = MAKE_TS(TA, TB, DA(L-1), DB(LEN-L-1)) ! construct TS
	CALL TEN_SAVE('TS',TEN_FLATTEN(TS,[1,3,5,7,0,2,4,6,8]))
	U = GS_SUBSPACE(TS, W)
!	TVAL = ANNEAL(TS, W) ! anneal W by TS
!	F = FA(L-1)+FB(LEN-L-1)+LOG(TVAL) ! calculate F
END SUBROUTINE LDMRG
! estimate iDMG trial state
FUNCTION MAKE_W5(WA, WB, P, PI) RESULT(W)
! call by IDMRG
	TYPE(TENSOR), INTENT(IN) :: WA, WB, P, PI
	TYPE(TENSOR) :: W
	
	! W = 
	!               2                 4
	!               │                 │
	!               2                 2
	! 1 ─ 1 P 2 ─ 3 WB 1 ─ 1 PI 2 ─ 1 WA 3 ─ 2 P 1 ─ 3
	W = TEN_PROD(TEN_PROD(TEN_PROD(P,WB,[2],[3]),PI,[2],[1]),TEN_PROD(P,WA,[2],[3]),[3],[2])
END FUNCTION MAKE_W5
! estimate fDMRG trial state
FUNCTION MAKE_W3(P, W1, W2) RESULT(W)
! called by IDMRG
! W1, W2 tensors can be either A or B type
	TYPE(TENSOR), INTENT(IN) :: P, W1, W2
	TYPE(TENSOR) :: W
	
	! W = 
	!               2        4
	!               │        │
	!               2        2
	! 1 ─ 1 P 2 ─ 3 W1 1 ─ 3 W2 1 ─ 3
	W = TEN_PROD(TEN_PROD(P,W1,[2],[3]),W2,[2],[3])
END FUNCTION MAKE_W3
! construct system tensor
FUNCTION MAKE_TS(TA, TB, DA, DB) RESULT (TS)
! input: TA,TB - site MPO, DA - A-block MPO, DB - B-block MPO
! output: TS - system MPO
! DB is optional: present => bulk, missing => boundary 
	TYPE(TENSOR), INTENT(IN) :: TA, TB, DA
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: DB
	TYPE(TENSOR) :: TS
	
	IF (PRESENT(DB)) THEN ! bulk algorithm
		! TS (bulk) = 
		!     1        3        7        5
		!     │        │        │        │
		!     3        3        3        3
		! ╭ 2 DA 1 ─ 2 TA 1 ─ 1 TB 2 ─ 1 DB 2 ╮
		! │   4        4        4        4    │
		! │   │        │        │        │    │
		! │   2        4        8        6    │
		! ╰───────────────────────────────────╯
		TS = TEN_PROD(TEN_PROD(DA,TA,[1],[2]),TEN_PROD(DB,TB,[1],[2]),[1,4],[1,4])
	ELSE ! boundary algorithm
		! TS (boundary) = 
		! 1        3        7             5
		! │        │        │             │
		! 3        3        3             1
		! DA 1 ─ 2 TA 1 ─ 1 TB 2 ─ 2 DA ⊗ I
		! 4        4        4             2
		! │        │        │             │
		! 2        4        8             6
		TS = TEN_PROD(TEN_PROD(TEN_PROD(DA,TA,[1],[2]),EYE_TEN([1,1])),TB,[1,4],[2,1])
	END IF
END FUNCTION MAKE_TS
! update DX given TX and WX
FUNCTION NEW_DX(TX, WX, DX) RESULT (DX1)
! input: TX - MPO, WX - MPS, DX - block MPO
! output: DX1 - new block MPO by packing TX and WX into the old
! DX optional: if present => bulk, missing => boundary
! zipper-order contraction algorithm is used
	TYPE(TENSOR), INTENT(IN) :: TX, WX
	TYPE(TENSOR), OPTIONAL, INTENT(IN) :: DX
	TYPE(TENSOR) :: DX1
	
	IF (PRESENT(DX)) THEN ! bulk algorithm
		! DX1 (bulk) = 
		!       ╭──── 1 WX* 3 ─ 3
		!       │        2
		!       │        │
		!       3        3
		! 2 ─ 2 DX 1 ─ 2 TX 1 ─ 1
		!       4        4
		!       │        │
		!       │        2 
		!       ╰───── 1 WX 3 ─ 4
		DX1 = TEN_PROD(TX,TEN_PROD(TEN_PROD(DX,TEN_CONJG(WX),[3],[1]),WX,[3],[1]),[2,3,4],[1,3,5])
	ELSE ! boundary algorithm
		! DX1 (boundary) = 
		!   ╭ 1 WX* 3 ─ 3
		!   │    2
		!   │    │
		!   ╯    3
		! 2 ── 2 TX 1 ─ 1
		!   ╮    4
		!   │    │
		!   │    2 
		!   ╰─ 1 WX 3 ─ 4
		DX1 = TEN_FLATTEN(TEN_PROD(TEN_PROD(TX,TEN_CONJG(WX),[3],[2]),WX,[3],[2]), [1,0,2,3,5,0,4,0,6])
	END IF
END FUNCTION NEW_DX
! rescale DX
SUBROUTINE RESCALE(DX, FX0, FX)
! input: DX - unscaled block MPO, FX0 - last level
! output: DX - scaled block MPO, FX - this scaling level
	USE CONST
	TYPE(TENSOR), INTENT(INOUT) :: DX
	COMPLEX, INTENT(IN)  :: FX0
	COMPLEX, INTENT(OUT) :: FX
	! local variables
	COMPLEX :: Z
	
	! get the scale of D by evaluation
	Z = EVALUATE(DX)
	IF (ABS(Z) /= 0.) THEN
		! rescale by a factor of Z
		DX%VALS = DX%VALS/Z
		! use it to calculate new level
		FX = FX0 + LOG(Z)
	END IF
END SUBROUTINE RESCALE
! evaluate environment tensor
FUNCTION EVALUATE(E) RESULT(Z)
! input: E - environment tensor
! output: Z - measurement value of E, by tensor trace
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: E
	COMPLEX :: Z
	! local tensor
	TYPE(TENSOR) :: E0
	
	! trace out the legs of E
	!  ╭───╮
	!  │   3
	! ╭┼ 1 E 2 ╮
	! ││   4   │
	! │╰───╯   │
	! ╰────────╯
	E0 = TEN_TRACE(E,[1,3],[2,4])
	! check that E0 is a number
	IF (SIZE(E0%DIMS) /= 0) THEN
		WRITE (*,'(A)') 'EVALUATE::xleg: not evaluated to a number.'
		STOP
	END IF
	! prepare output
	IF (SIZE(E0%VALS) == 0) THEN ! if E0 empty
		Z = Z0 ! return zero
	ELSE ! if not empty
		Z = E0%VALS(1) ! return the value
	END IF
END FUNCTION EVALUATE
! cal entanglement entropy
FUNCTION ENTROPY(P) RESULT (S)
! input: P - Schmidt spectrum
! update: S - entanglement entropy
	TYPE(TENSOR), INTENT(IN) :: P
	REAL :: S
	! local variables
	REAL, ALLOCATABLE :: EVALS(:)
	
	EVALS = REALPART(P%VALS)**2
	S = -SUM(EVALS*LOG(EVALS))
END FUNCTION ENTROPY
! ----------- Solvers ------------
! anneal the state W to the ground state of TS
FUNCTION ANNEAL(TS, W) RESULT (TVAL)
! input: TS - system transfer tensor, W - state
! on output: W  is modified to the fixed point state
! return the corresponding eigen value TVAL
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	
! +++++++++ choose a solver here +++++++++++
!	TVAL = ANNEAL0(TS, W) ! self-made (fas, unstable)
	TVAL = ANNEAL1(TS, W) ! package (slower, robust)
! ++++++++++++++++++++++++++++++++++++++++++
END FUNCTION ANNEAL
! self-made Arnoldi algorithm (fast, unstable)
FUNCTION ANNEAL0(TS, W) RESULT (TVAL)
! called by ANNEAL
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
	INTEGER, PARAMETER :: N = 16 ! Krylov space dimension
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	! local variables
	INTEGER :: DIM, I, J, K, ITER, INFO
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:)
	COMPLEX, ALLOCATABLE :: V(:,:)
	COMPLEX :: TVAL0, A(N,N), D(N),Z(N,N),WORK(65*N)
	REAL :: RWORK(2*N)
	LOGICAL :: EXH
	
	! unpack data from tensor
	! collect leg-combined inds in TS and W (remember to +1)
	LINDS = COLLECT_INDS(TS,[1,3,5,7])+1
	RINDS = COLLECT_INDS(TS,[2,4,6,8])+1
	WINDS = COLLECT_INDS(W,[1,2,3,4])+1
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(V(DIM,N))
	V = Z0 ! clear
	! dump W to the 1st col of V
	FORALL (I = 1:SIZE(W%INDS))
		V(WINDS(I),1) = W%VALS(I)
	END FORALL
	V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
	! prepare to start Arnoldi iteration
	TVAL0 = Z0
	EXH = .FALSE. ! space exhausted flag
	! use Arnoldi iteration algorithm
	DO ITER = 1, MAX_ITER
		A = Z0 ! initialize Heisenberg matrix
		! construct Krylov space
		DO K = 2, N
			! apply TS to V(:,K-1) -> V(:,K)
			V(:,K) = Z0
			DO I = 1,SIZE(TS%INDS)
				V(LINDS(I),K) = V(LINDS(I),K) + TS%VALS(I)*V(RINDS(I),K-1)
			END DO
			! orthogonalization by stabilized Gram–Schmidt process
			DO J = 1, K-1
				A(J,K-1) = DOT_PRODUCT(V(:,J),V(:,K))
				V(:,K) = V(:,K) - A(J,K-1)*V(:,J)
			END DO
			! cal the norm of residual vector
			A(K,K-1) = SQRT(DOT_PRODUCT(V(:,K),V(:,K)))
			! if it is vanishing, the Arnoldi iteration has broken down
			IF (ABS(A(K,K-1))<TOL) THEN
				EXH = .TRUE.
				EXIT ! exit the iteration, stop construction of Krylov space
			END IF
			! otherwise, normalize the residual vect to a new basis vect
			V(:,K) = V(:,K)/A(K,K-1)
		END DO !K
		! now the Heisenberg matrix has been constructed
		! the action of TS is represented on the basis V as A
		! call LAPACK to diagonalize A
		CALL ZGEEV('N','V',N,A,N,D,Z,1,Z,N,WORK,65*N,RWORK,INFO)
		! now D holds the eigen vals, and Z the eigen vects
		! find the max abs eigen val
		I = LOC_MAX_MODE(D)
		TVAL = D(I) ! save it in Tval
		! reorganize the eigen vector from Z(:,I)
		V(:,1) = MATMUL(V,Z(:,I)) ! save to 1st col of V
		V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
		! check convergence
		! if space exhausted, or relative error < tol
		IF (EXH .OR. ABS((TVAL-TVAL0)/TVAL) < TOL) THEN
			! Arnoldi iteration has converge
			EXIT ! exit Arnoldi interation
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
		END IF
	END DO ! next Arnoldi interation
	! if exceed max iteration
	IF (ITER > MAX_ITER) THEN !then power iteration has not converge
		WRITE (*,'(A)') 'ANNEAL::fcov: Arnoldi iteration failed to converge.'
	END IF
	! reconstruct W tensor for output
	W%INDS = [(I,I=0,DIM-1)]
	W%VALS = V(:,1)
END FUNCTION ANNEAL0
! by calling to LAPACK and ARPACK (slower, robust)
FUNCTION ANNEAL1(TS, W) RESULT (TVAL)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS
	TYPE(TENSOR), INTENT(INOUT) :: W
	COMPLEX :: TVAL
	! parameters
    CHARACTER(*),PARAMETER :: BMAT = 'I'    ! regular eigen problem
    CHARACTER(*),PARAMETER :: WHICH = 'LM'  ! largest magnitude mode
    CHARACTER(*),PARAMETER :: HOWMNY = 'A'  ! get Ritz eigen vecs
    LOGICAL, PARAMETER     :: RVEC = .TRUE. ! calculate eigen vecs
    INTEGER, PARAMETER :: MAXITER = 300 ! max Arnoldi iterations
    REAL, PARAMETER    :: TOL = 1.D-12 ! tolerance
    INTEGER, PARAMETER :: NCV = 16 ! Krylov vector space dim
    INTEGER, PARAMETER :: NEV = 6  ! num of eigenvalues to be compute
    INTEGER, PARAMETER :: LWORKL = 3*NCV**2 + 5*NCV
	! local variables
	INTEGER :: NCALL ! num of calls of TS*W by ARPACK
	CHARACTER, PARAMETER :: BACK(20) = CHAR(8)
	! tensor related
	INTEGER :: I, P ! iterator/index
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:) ! inds to collect
	INTEGER :: N ! Hilbert space dim
	COMPLEX, ALLOCATABLE :: D(:)   ! eigen values
	COMPLEX, ALLOCATABLE :: Z(:,:) ! eigen vectors
	COMPLEX, ALLOCATABLE :: A(:,:) ! mat rep of TS
	! ARPACK related
	INTEGER :: IDO, INFO ! communication/information flag
	INTEGER :: NCONV  ! num of converged eigen vals
	INTEGER :: IPARAM(11), IPNTR(14) ! used by ARPACK
	COMPLEX :: SIGMA  ! shift, not referenced
	LOGICAL, ALLOCATABLE :: SEL(:) ! eigen vect selection
	REAL,    ALLOCATABLE :: RWORK(:) ! work space
	COMPLEX, ALLOCATABLE :: RESID(:) ! residual vect (to put initial vect)
	COMPLEX, ALLOCATABLE :: V(:,:)   ! Arnoldi basis vectors
	COMPLEX, ALLOCATABLE :: WORKL(:), WORKEV(:) ! work spaces
	COMPLEX, ALLOCATABLE, TARGET :: WORKD(:) ! work space
	COMPLEX, POINTER :: W0(:), W1(:) ! pointer to WORKD
	
	! cal total dim of W
	N = PRODUCT(W%DIMS) ! this will be the Hilbert space dim
	IF (N <= NCV) THEN ! for small scale problem
		! call LAPACK directly
		! unpack data from TS
		A = TEN2MAT(TS,[1,3,5,7],[2,4,6,8])
		IF (.NOT.(SIZE(A,1) == N .AND. SIZE(A,2) == N)) THEN
			WRITE (*,'(A)') 'ANNEAL::xdim: incompatible dimension between TS and W.'
			STOP
		END IF
		! allocate for work space
		ALLOCATE(D(N), Z(N,N), WORKD(65*N), RWORK(2*N))
		! call LAPACK to diagonalize A
		CALL ZGEEV('N','V', N, A, N, D, Z, 1, Z, N, WORKD,&
			65*N, RWORK, INFO)
		IF (INFO /= 0) THEN ! error handling
			IF (INFO > 0) THEN
				WRITE (*,'(A)') 'ANNEAL::fcnv: the QR algorithm failed.'
			ELSE
				WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZGEEV error ', INFO
			END IF
		END IF
		! now D holds the eigen vals, and Z the eigen vects
		! find the max abs eigen val
		P = LOC_MAX_MODE(D)
	ELSE ! for large scale problem
		! call ARPACK
		! allocate for work space
		ALLOCATE(D(NCV), Z(N,NEV), RESID(N), V(N,NCV), WORKD(3*N),&
			WORKL(LWORKL), WORKEV(3*NCV), RWORK(NCV), SEL(NCV), STAT = INFO)
		! collect leg-combined inds in TS and W (remember to +1)
		LINDS = COLLECT_INDS(TS,[1,3,5,7])+1
		RINDS = COLLECT_INDS(TS,[2,4,6,8])+1
		WINDS = COLLECT_INDS(W,[1,2,3,4])+1
		! unpack initial vect from W tensor input
		RESID = Z0 ! clear
		FORALL (I = 1:SIZE(W%INDS))
			RESID(WINDS(I)) = W%VALS(I) ! dump
		END FORALL
		RESID = RESID/SQRT(DOT_PRODUCT(RESID,RESID)) ! normalize
		! set ARPACK parameters
		IPARAM(1) = 1       ! ishift
		IPARAM(3) = MAXITER ! maxitr
		IPARAM(7) = 1       ! model 1: A*x = lambda*x
		! now it is ready to start Arnoldi iteration
		IDO  = 0 ! initialize communication flag
		INFO = 1 ! RESID contains the initial residual vector
		! call ARPACK ----------------------
		! main loop (reverse communication)
		NCALL = 0
		WRITE (6,'(A,I5,A)',ADVANCE = 'NO') 'dim:',N,', iter:'
		DO WHILE (.TRUE.) ! repeatedly call ZNAUPD
			! until converge or exceed MAXITER
			CALL ZNAUPD(IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, &
				V, N, IPARAM, IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO)
			! take action according to communication flag IDO
			SELECT CASE (IDO)
				CASE (-1,1) ! compute W1 = TS*W0
					! monitor each call
					NCALL = NCALL +1 ! NCALL inc
					WRITE (6,'(I4)',ADVANCE = 'NO') NCALL
					! perform action of TS to W0 -> W1
					! point W0, W1 to the correct part of WORKD
					W0 => WORKD(IPNTR(1):IPNTR(1)+N-1)
					W1 => WORKD(IPNTR(2):IPNTR(2)+N-1)
					! carry out sparse matrix multiplication
					W1 = Z0 ! clear
					DO I = 1,SIZE(TS%INDS)
						W1(LINDS(I)) = W1(LINDS(I)) + TS%VALS(I)*W0(RINDS(I))
					END DO
					! monitor back space
					WRITE (6,'(4A1)',ADVANCE = 'NO') BACK(1:4)
				CASE (99) ! done
					EXIT ! exit Arnoldi iteration
			END SELECT
		END DO
		! clear the monitor
		WRITE (6,'(20A1)',ADVANCE = 'NO') BACK
		! Either we have convergence or there is an error    
		IF (INFO /= 0) THEN ! error handling
			SELECT CASE (INFO)
				CASE (1)
					WRITE (*,'(A)') 'ANNEAL::fcnv: Arnoldi iteration failed to converge.'
				CASE (3)
					WRITE (*,'(A)') 'ANNEAL::xsft: no shift can be applied, try larger NCV.'        	
				CASE DEFAULT
					WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZNAUPD error ', INFO
					STOP
			END SELECT
		ELSE
			! No fatal errors occurred.
			! Post-Process using ZNEUPD.
			CALL ZNEUPD(RVEC, HOWMNY, SEL, D, Z, N, SIGMA, WORKEV, &
				BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, N, IPARAM, &
				IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO)
			IF (INFO /=0 ) THEN ! error handling
				WRITE (*,'(A,I3)') 'ANNEAL::xinf: ZNEUPD error ', INFO
				STOP
			END IF
			NCONV = IPARAM(5) ! get the num of convergent modes
		END IF
		! now D(:NCONV) holds the eigen vals
		! and Z(:,:NCONV) the eigen vects
		! find the max abs eigen val
		P = LOC_MAX_MODE(D(1:NCONV)) ! get position
	END IF
	! prepare output
	TVAL = D(P) ! output eigen val as TVAL
	! reconstruct W tensor from the corresponding mode
	W%INDS = [(I,I=0,N-1)] ! make index
	W%VALS = Z(:,P) ! dump data
END FUNCTION ANNEAL1
! locate max mode from eigen values
FUNCTION LOC_MAX_MODE(D) RESULT (P)
! input: D - eigen values of TS
! output: P - position where the max mode is located
	COMPLEX, INTENT(IN) :: D(:)
	INTEGER :: P
	! local variables
	INTEGER :: I
	COMPLEX :: LGDP, LGDI
	COMPLEX, ALLOCATABLE :: LGD(:)
	REAL, PARAMETER :: TOL = 1.E-7 ! allowed relative error
	
	! take log of D
	LGD = LOG(D)
	! assuming  P at the beginning
	P = 1
	LGDP = LGD(P) ! record the in-hand lgD
	! start searching for better choice
	DO I = 2, SIZE(D)
		LGDI = LGD(I) ! get the current lgD
		! compare with the one in hand
		IF (REAL(LGDI) > REAL(LGDP)+TOL) THEN
			! lgD(I) > lgD(P)
			P = I ! catch I
			LGDP = LGD(P) ! record in-hand
		ELSEIF (REAL(LGDI) > REAL(LGDP)-TOL) THEN
			! lgD(I) ~ lgD(P) (degenerated in magnitude)
			IF (ABS(IMAG(LGDI)) < ABS(IMAG(LGDP)) - TOL) THEN
				! |phase of I| < |phase of P|, D(I) is more real
				P = I ! catch I
				LGDP = LGD(P) ! record in-hand
			ELSEIF (ABS(IMAG(LGDI)) < ABS(IMAG(LGDP))+TOL) THEN
				! |phase of I| ~ |phase of P| (phase degenerated)
				IF (IMAG(LGDI) > IMAG(LGDP)) THEN
					! choose the positive one
					P = I ! catch I
					LGDP = LGD(P) ! record in-hand
				END IF
			END IF
		END IF
	END DO
END FUNCTION LOC_MAX_MODE
! find ground state subspace
FUNCTION GS_SUBSPACE(TS, W) RESULT (U)
	USE CONST
	TYPE(TENSOR), INTENT(IN) :: TS, W
	TYPE(TENSOR) :: U
	! parameters
	INTEGER, PARAMETER :: N = 16 ! Krylov space dimension
	INTEGER, PARAMETER :: MAX_ITER = 500 ! max interation
	REAL, PARAMETER :: TOL = 1.E-12 ! allowed error of Tval
	! local variables
	INTEGER :: DIM, I, J, K, ITER, INFO
	INTEGER, ALLOCATABLE :: LINDS(:), RINDS(:), WINDS(:)
	COMPLEX, ALLOCATABLE :: V(:,:)
	COMPLEX :: TVAL, TVAL0, A(N,N), D(N),Z(N,N),WORK(65*N)
	REAL :: RWORK(2*N)
	LOGICAL :: EXH
	
	! unpack data from tensor
	! collect leg-combined inds in TS and W (remember to +1)
	LINDS = COLLECT_INDS(TS,[1,3,5,7])+1
	RINDS = COLLECT_INDS(TS,[2,4,6,8])+1
	WINDS = COLLECT_INDS(W,[1,2,3,4])+1
	! cal total dim of W
	DIM = PRODUCT(W%DIMS)
	! allocate Krylov space
	ALLOCATE(V(DIM,N))
	V = Z0 ! clear
	! dump W to the 1st col of V
	FORALL (I = 1:SIZE(W%INDS))
		V(WINDS(I),1) = W%VALS(I)
	END FORALL
	V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
	! prepare to start Arnoldi iteration
	TVAL0 = Z0
	EXH = .FALSE. ! space exhausted flag
	! use Arnoldi iteration algorithm
	DO ITER = 1, MAX_ITER
		A = Z0 ! initialize Heisenberg matrix
		! construct Krylov space
		DO K = 2, N
			! apply TS to V(:,K-1) -> V(:,K)
			V(:,K) = Z0
			DO I = 1,SIZE(TS%INDS)
				V(LINDS(I),K) = V(LINDS(I),K) + TS%VALS(I)*V(RINDS(I),K-1)
			END DO
			! orthogonalization by stabilized Gram–Schmidt process
			DO J = 1, K-1
				A(J,K-1) = DOT_PRODUCT(V(:,J),V(:,K))
				V(:,K) = V(:,K) - A(J,K-1)*V(:,J)
			END DO
			! cal the norm of residual vector
			A(K,K-1) = SQRT(DOT_PRODUCT(V(:,K),V(:,K)))
			! if it is vanishing, the Arnoldi iteration has broken down
			IF (ABS(A(K,K-1))<TOL) THEN
				EXH = .TRUE.
				EXIT ! exit the iteration, stop construction of Krylov space
			END IF
			! otherwise, normalize the residual vect to a new basis vect
			V(:,K) = V(:,K)/A(K,K-1)
		END DO !K
		! now the Heisenberg matrix has been constructed
		! the action of TS is represented on the basis V as A
		! call LAPACK to diagonalize A
		CALL ZGEEV('N','V',N,A,N,D,Z,1,Z,N,WORK,65*N,RWORK,INFO)
		! now D holds the eigen vals, and Z the eigen vects
		! find the max abs eigen val
		I = LOC_MAX_MODE(D)
		TVAL = D(I) ! save it in Tval
		! reorganize the eigen vector from Z(:,I)
		V(:,1) = MATMUL(V,Z(:,I)) ! save to 1st col of V
		V(:,1) = V(:,1)/SQRT(DOT_PRODUCT(V(:,1),V(:,1))) ! normalize
		! check convergence
		! if space exhausted, or relative error < tol
		IF (EXH .OR. ABS((TVAL-TVAL0)/TVAL) < TOL) THEN
			! Arnoldi iteration has converge
			EXIT ! exit Arnoldi interation
		ELSE ! not converge, next iteration
			TVAL0 = TVAL ! save TVAL to TVAL0
		END IF
	END DO ! next Arnoldi interation
	! if exceed max iteration
	IF (ITER > MAX_ITER) THEN !then power iteration has not converge
		WRITE (*,'(A)') 'ANNEAL::fcov: Arnoldi iteration failed to converge.'
	END IF
	! reconstruct W tensor for output
!	W%INDS = [(I,I=0,DIM-1)]
!	W%VALS = V(:,1)
END FUNCTION GS_SUBSPACE
! ---------- grow column tensor -------------
! grow column tensor
FUNCTION GROW(T, L) RESULT (M)
! input: T - MPO tensor, L - layers to grow
! output: M - L layers of MPO tensors
! staking 1 M 2 =  ─ 1 T 2 ─ 1 T 2 ─ ...
	TYPE(TENSOR), INTENT(IN) :: T
	INTEGER :: L
	TYPE(HUGE_TENSOR) :: M
	! local variables
	INTEGER :: I
	
	M = TOHUGE(T) ! put the first layer down
	DO I = 2, L
		!           3
		!       ╭───┴───╮
		!       3       3
		! 1 ─ 1 M 2 ─ 1 T 2 ─ 2
		!       4       4
		!       ╰───┬───╯
		!           4
		M%TEN = SVD_REDUCE(TEN_FLATTEN(TEN_PROD(M%TEN,T,[2],[1]),[1,0,4,0,2,5,0,3,6]))
		CALL NORMALIZE(M)
		WRITE (*,'(I3,A,F5.1,A)') I,':',SVD_ERR*100,'%'
!		WRITE (*,'(F8.4)') M%LEV/LOG(2.)
	END DO
END FUNCTION GROW
! reduce tensor by SVD
FUNCTION SVD_REDUCE(T1) RESULT (T2)
! input: T1 to be reduct, output: T2 as the result
! reduce leg 1,2 by HOSVD
	USE MODEL
	TYPE(TENSOR), INTENT(IN) :: T1
	TYPE(TENSOR) :: T2
	! local tensors
	TYPE(TENSOR), ALLOCATABLE :: US(:)
	INTEGER, PARAMETER :: DCUT = 128
	
	T2 = T1 ! copy T1 to T2
	! perform HOSVD
	CALL HOSVD(T2, US, [3,4], [DCUT,DCUT])
	!       │
	!       3
	! 1 ─ 1 T2 2 ─ 2
	!       4
	!       ╰─ 2 US2 1 ━ 1 US1 2 ─ 4
	T2 = TEN_PROD(T2,TEN_PROD(US(2),US(1),[1],[1]),[4],[1])
END FUNCTION SVD_REDUCE
! ----------- measurement -------------------
! ----------- Debug ------------
! show status after DMRG step
SUBROUTINE SHOW_DMRG_STATUS(MODE,L,F,P)
	USE MODEL
	CHARACTER, INTENT(IN) :: MODE
	INTEGER, INTENT(IN) :: L
	COMPLEX, INTENT(IN) :: F
	TYPE(TENSOR), INTENT(IN) :: P
	! local variables
	REAL :: S
	INTEGER :: I1, I2
	CHARACTER(2) :: JN
	
	SELECT CASE (MODE)
		CASE ('I')
			I1 = L
			I2 = L+1
			JN = '<>'
		CASE ('F')
			I1 = L
			I2 = L+1
			JN = '->'
			IF (I2 == LEN+1) I2 = 1
		CASE ('B')
			I1 = LEN-L
			I2 = LEN-L+1
			JN = '<-'
			IF (I1 == 0) I1 = LEN
	END SELECT
	IF (L == 0 .OR. L == LEN) THEN
		S = 0.
	ELSE
		S = ENTROPY(P)
	END IF
	WRITE (*,'(I3,A,I3,A,F10.6,A,F6.3,A,F10.6,A,F5.2,A)') I1,JN,I2,': F = (', REALPART(F),',', IMAGPART(F), '), S = ', S/LOG(2.), 'bit ', SVD_ERR*100,'%'
END SUBROUTINE SHOW_DMRG_STATUS
! reconstruct wave function
FUNCTION MPS(WA) RESULT (W)
! input: WA - MPS tensors
! output: W - MPS state vector
	TYPE(TENSOR), INTENT(IN) :: WA(:)
	COMPLEX, ALLOCATABLE :: W(:)
	! local variables
	TYPE(TENSOR) :: WT
	INTEGER :: NW, IW, I
	
	NW = SIZE(WA) ! get size of WA
	IF (NW == 0) THEN
		W = [COMPLEX::]
	END IF
	WT = WA(1)
	DO IW = 2, NW
		WT = TEN_FLATTEN(TEN_PROD(WT, WA(IW),[3],[1]),[1,0,2,3,0,4])
	END DO
	WT = TEN_TRACE(WT,[1],[3])
	ALLOCATE(W(WT%DIMS(1)))
	FORALL (I = 1:SIZE(WT%INDS))
		W(WT%INDS(I)+1) = WT%VALS(I)
	END FORALL
END FUNCTION MPS
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
CONTAINS
! ------------ Data --------------
! collect data
SUBROUTINE COLLECT(BETAS)
	USE MODEL
	REAL, INTENT(IN) :: BETAS(:) ! a list of beta to be sampled
	! local tensors
	TYPE(TENSOR) :: T, WS(LEN)
	TYPE(TENSOR) :: M, OS(2)
	! local variables
	INTEGER :: I, N
	
	! prepare observables to be measured
!	DO I = 1,2 ! two operators are both sigma_3
!		OS(I) = PAULI_MAT([3])
!	END DO
!	N = SIZE(BETAS) ! get size of beta list
!	DO I = 1, N ! for each beta in the list
!		BETA = BETAS(I) ! (remember to set beta!!!)
!		CALL SET_MPO(T) ! update MPO, using new beta
!		! prepare to launch DMRG
!		WRITE (*,'(A,I3,A,F5.2,A,F5.2,A,F5.2)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA, ', crossing = ', CROSS
!		! launch DMRG to find ground state MPS given MPO
!		CALL DMRG(T, WS)
!		! take measurements on the MPS
!		M = MEASURE(WS, OS)
!		! save result to disk
!		CALL TEN_SAVE('./data center (IS)/M'//FILENAME(),M)
!	END DO 
END SUBROUTINE COLLECT
! make filename
FUNCTION FILENAME()
! file naming system
! call by COLLECT, will use MODEL parameters
	USE MODEL
	CHARACTER(17) :: FILENAME
	
	WRITE (FILENAME,'(A,I3.3,A,I1,A,I1,A,I1,A,SP,I5.4)') 'D',MAX_CUT,'S',SWEEPS,'Q',NINT(THETA/PI),'C',NINT(CROSS),'B',NINT(BETA*10000)
END FUNCTION FILENAME
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	INTEGER :: N, I
	
	I = 131072
	N = I*I
	PRINT *, N
END SUBROUTINE TEST
! test DMRG
SUBROUTINE TEST_DMRG()
	USE MATHIO
	USE MODEL
	TYPE(TENSOR) :: TB
	TYPE(HUGE_TENSOR) :: M
	REAL :: T0, T1
	
	WRITE (*,'(A,I3,A,F5.2,A,F5.2,A,F5.2)') 'cut = ', MAX_CUT, ', theta = ', THETA/PI, '*pi, beta = ', BETA, ', crossing = ', CROSS
	CALL CPU_TIME(T0)
	CALL SET_MPO(TB)
	CALL DMRG(TB, M)
	CALL CPU_TIME(T1)
!	PRINT *, T1-T0
!	M%TEN%VALS = M%TEN%VALS * EXP(M%LEV)
!	CALL TEN_SAVE('M', M%TEN)
END SUBROUTINE TEST_DMRG
! test transfer matrix method
SUBROUTINE TEST_TRANSF()
	USE MATHIO
	INTEGER, PARAMETER :: L = 4
	TYPE(TENSOR) :: T, A
	INTEGER :: I
	COMPLEX, ALLOCATABLE :: M0(:,:), M3(:,:)
	
	CALL SET_MPO(T) ! obtain MPO tensor
	A = T ! set first site
	DO I = 2, L ! grow to length L
		!           3
		!       ╭───┴───╮
		!       3       3
		! 1 ─ 1 A 2 ─ 1 T 2 ─ 2
		!       4       4
		!       ╰───┬───╯
		!           4
		A = TEN_FLATTEN(TEN_PROD(A,T,[2],[1]),[1,0,4,0,2,5,0,3,6])
	END DO
	! put matrix
	M0 = TEN2MAT(TEN_TRACE(A,[1],[2]),[1],[2])
	M3 = TEN2MAT(TEN_PROD(A,PAULI_MAT([3]),[1,2],[2,1]),[1],[2])
	CALL EXPORT('M0',M0)
	CALL EXPORT('M3',M3)
END SUBROUTINE TEST_TRANSF
! test transfer matrix method
SUBROUTINE TEST_TRANSF2()
	USE MODEL
	USE MATHIO
	TYPE(TENSOR) :: T
	TYPE(HUGE_TENSOR) :: M
	INTEGER :: I
	COMPLEX, ALLOCATABLE :: M0(:,:), M3(:,:)
	
	CALL SET_MPO(T) ! obtain MPO tensor
	M = GROW(T, LEN)
	! put matrix
	M0 = TEN2MAT(TEN_TRACE(M%TEN,[1],[2]),[1],[2])*EXP(M%LEV)
	M3 = TEN2MAT(TEN_PROD(M%TEN,PAULI_MAT([3]),[1,2],[2,1]),[1],[2])*EXP(M%LEV)
	CALL EXPORT('M0',M0)
	CALL EXPORT('M3',M3)
END SUBROUTINE TEST_TRANSF2
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	INTEGER :: I
	PRINT *, '------------ DMRG -------------'

!	CALL TEST()
	CALL TEST_DMRG()
!	CALL TEST_TRANSF()
!	CALL TEST_TRANSF2()
!	CALL COLLECT([0.4406])
END PROGRAM MAIN