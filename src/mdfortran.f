      PROGRAM MDFORTRAN
      ! Program written by Tridip Das (dastridi)
         INTEGER :: TOT_STEPS, STEPS, L
         REAL :: DEL_TM, LAT, TGOAL, SIGMA, EPS, M, TIME
        TGOAL = 1.0
        SIGMA = 1.0
        EPS = 1.0
        M  = 1.0
        L = 6
        TOT_STEPS = 2000
        DEL_TM = 0.001
        LAT = 1.5874
      !
        CALL CPU_TIME(TIME)
        PRINT *, " Start time ", TIME
        CALL MDMAIN(TGOAL,SIGMA,EPS,M,L,LAT,TOT_STEPS,DEL_TM)
        CALL CPU_TIME(TIME)
        PRINT *, " End time ", TIME
      !
      END PROGRAM MDFORTRAN
      !
      ! The below subroutine was supposed to be called by md.   
      !
      SUBROUTINE MDMAIN(TGOAL,SIGMA,EPS,MASS,L,LAT,TOT_STEPS,DEL_TM)
      ! This is the main logic for MD simulation run
        IMPLICIT NONE
      ! Variable Declaration
        INTEGER :: N, I, J, K, TOT_STEPS, STEPS, L
        REAL :: DEL_TM, LAT, PI, VOLUME, DENSITY, VIRIAL, PRESSURE
        REAL, ALLOCATABLE :: ATOM_POS(:,:), VEL(:,:), ACCN(:,:)
      ! Individual atomic positions, Velocity, Acceleration
        REAL, ALLOCATABLE :: INDATA(:), E_POT(:), E_KIN(:), PRESS(:)
        REAL, ALLOCATABLE :: TEMPERATURE(:), MOMENTUM(:)
        REAL, ALLOCATABLE :: E_POT_AVG(:), E_KIN_AVG(:), E_TOT(:)
      ! Potential, Kinetic, & Total Energy
        REAL :: RIJ(3), RSQUARE, EPS, SIGMA
        REAL :: POT_EN, F(3), MASS, RUIJ(3),VEL_TOT(3),VEL_VEC(3),TGOAL
        REAL :: RAND1, RAND2, kB = 1.0, VUIJ(3), VEL_WRITE(5)
      ! 
        N = 4*L**3
        VOLUME = (LAT*L)**3
        DENSITY = N/VOLUME
      ! 
      ! Allocate variables
        ALLOCATE(ATOM_POS(N,3))
        ALLOCATE(VEL(N,3))
        ALLOCATE(ACCN(N,3))
        ALLOCATE(INDATA(N*3))
        ALLOCATE(E_POT(N))
        ALLOCATE(E_KIN(N))
        ALLOCATE(E_POT_AVG(TOT_STEPS))
        ALLOCATE(E_KIN_AVG(TOT_STEPS))
        ALLOCATE(E_TOT(TOT_STEPS))
        ALLOCATE(TEMPERATURE(TOT_STEPS))
        ALLOCATE(PRESS(TOT_STEPS))
        ALLOCATE(MOMENTUM(TOT_STEPS))
      !
      ! Read the initial atom positions
        OPEN(UNIT=11,FILE='initial_lattice.dat',STATUS='OLD',           &
     &       ACTION='READ')
        READ(11,*) INDATA
      !
      ! Open output file to write Energy related data
          OPEN(UNIT=31,FILE='final_lat.dat',STATUS='REPLACE',           &
     &       ACTION='WRITE')
      !
          OPEN(UNIT=32,FILE='Energy.dat',STATUS='REPLACE',              &
     &       ACTION='WRITE')
        WRITE(32,*) " Time ", " Potential_Energy ", " Kinetic_Energy ", &
     &              " Total_Energy ", " Temperature ", " Moment of Cen" &
     &              , " Pressure "
      !
          OPEN(UNIT=33,FILE='Velocity.dat',STATUS='REPLACE',            &
     &       ACTION='WRITE')
      !
      ! Open output file to write final lattice positions 
 
      !  
        DO I = 1, N 
           DO J = 1, 3
              ATOM_POS(I,J) = INDATA((I-1)*3+J)          
           ENDDO
      !     print *, ATOM_POS(I,:)
        ENDDO  
      !
      ! Initialize velocity based on temperature and kinetic
      ! energy
        PI = 3.141592653
      ! Box Muller algorithm
        DO I = 1, N
           CALL RANDOM_NUMBER(RAND1)
           CALL RANDOM_NUMBER(RAND2)
           VEL(I,1) = SQRT(-2*LOG(RAND1))*COS(2*PI*RAND2)
           VEL(I,2) = SQRT(-2*LOG(RAND1))*SIN(2*PI*RAND2)
           CALL RANDOM_NUMBER(RAND1)
           CALL RANDOM_NUMBER(RAND2)
           VEL(I,3) = SQRT(-2*LOG(RAND1))*COS(2*PI*RAND2)
      !     VUIJ = VEL(I,:)/DOT_PRODUCT(VEL(I,:),VEL(I,:))
              VEL(I,:) = VEL(I,:)*sqrt(3*kB*TGOAL*(N-1)/N)
      !     print *, VEL(I,:)
        ENDDO
      ! Scale the velocity to keep centre of mass constant
        DO J = 1, 3
           VEL_TOT(J) = SUM(VEL(:,J))
        ENDDO
      !
        DO I = 1, N
           VEL(I,:) = VEL(I,:) - (VEL_TOT/N)
      !     print *, VEL(I,:)
        ENDDO  
      ! End of velocity Initialization
      !
      ! Main loop
      TIME: DO STEPS = 1, TOT_STEPS
      !
      ! Force calculation
        CALL FORCE(N,ATOM_POS,RIJ,L,LAT,E_POT,ACCN,VIRIAL) 
      !
      ! Calculate Kinetic Energy and Temperature   
           DO I = 1, N
              VEL_VEC = VEL(I,:)
              E_KIN(I) = 0.5*MASS*DOT_PRODUCT(VEL_VEC, VEL_VEC) ! Kinetic energy of each particle
           ENDDO
      !
           E_KIN_AVG(STEPS) = SUM(E_KIN)/N ! Average K.E. per particle
           TEMPERATURE(STEPS) = 2*E_KIN_AVG(STEPS)/3
      !
      !    Perform velocity scaling
           IF (MOD(STEPS,20) == 0) THEN
              DO I = 1, N
                 VEL(I,:) = VEL(I,:)*SQRT(TGOAL/TEMPERATURE(STEPS))
              ENDDO
           ENDIF
      !
           VEL = VEL + 0.5 * DEL_TM * ACCN
      !
           ATOM_POS = ATOM_POS +DEL_TM *VEL + 0.5*(DEL_TM**2)*ACCN
      !
           CALL FORCE(N,ATOM_POS,RIJ,L,LAT,E_POT,ACCN,VIRIAL)
      !
           VEL = VEL + 0.5 * DEL_TM * ACCN
      !
      ! update velocity for momentum calculation
           DO J = 1, 3
              VEL_TOT(J) = SUM(VEL(:,J))
           ENDDO
      !
           DO I = 1, N
              VEL(I,:) = VEL(I,:) - (VEL_TOT/N)
           ENDDO
      !
      ! Momentum of the centre of mass
           DO J = 1, 3
              VEL_TOT(J) =SUM(VEL(:,J))
           ENDDO
           MOMENTUM(STEPS) = MASS*SQRT(DOT_PRODUCT(VEL_TOT,VEL_TOT))
      !     
           ATOM_POS = ATOM_POS +DEL_TM *VEL + 0.5*(DEL_TM**2)*ACCN
      !    Apply periodic boundary condition
           DO I = 1, N
              DO J = 1, 3
                 IF (ATOM_POS(I,J) > L*LAT) THEN
                     ATOM_POS(I,J) = ATOM_POS(I,J) - L*LAT
                 ELSEIF (ATOM_POS(I,J) < 0.0) THEN
                     ATOM_POS(I,J) = ATOM_POS(I,J) + L*LAT
                 ELSE
                    CONTINUE
                 ENDIF
              ENDDO
           ENDDO
      !
           E_POT_AVG(STEPS) = SUM(E_POT)/N ! Avg Pot En per particles
      !
      ! Calculate Kinetic Energy and Temperature
           DO I = 1, N
              VEL_VEC = VEL(I,:)
              E_KIN(I) = 0.5*MASS*DOT_PRODUCT(VEL_VEC, VEL_VEC) ! Kinetic energy of each particle
           ENDDO
      !
           E_KIN_AVG(STEPS) = SUM(E_KIN)/N ! Average K.E. per particle
   
      ! Pressure of each step
           PRESS(STEPS) = DENSITY*TEMPERATURE(STEPS) + VIRIAL/VOLUME
      ! 
           E_TOT(STEPS) = E_POT_AVG(STEPS) + E_KIN_AVG(STEPS)
           WRITE(32,*) (STEPS * DEL_TM), E_POT_AVG(STEPS),              &
     &               E_KIN_AVG(STEPS), E_TOT(STEPS),TEMPERATURE(STEPS), &
     &               MOMENTUM(STEPS),PRESS(STEPS)
      ! 
      ! Write velocity of first five particles with time
           DO I = 1, 5
              VEL_WRITE(I) = SQRT(DOT_PRODUCT(VEL(I,:),VEL(I,:)))
           ENDDO
           WRITE(33,*) (STEPS * DEL_TM), VEL_WRITE
      !
      ENDDO TIME
      ! 
        DO I = 1, N
           WRITE (31,*) ATOM_POS(I,:)
        ENDDO
      !        
        CLOSE(11)
        CLOSE(31)
        CLOSE(32)
        CLOSE(33)
      ! 
      END SUBROUTINE MDMAIN
      !
      SUBROUTINE FORCE(N,ATOM_POS,RIJ,L,LAT,E_POT,ACCN,VIRIAL)
      !
        INTEGER :: N, I, J, K, L
        REAL :: ATOM_POS(N,3), ACCN(N,3), E_POT(N)
      ! Individual atomic positions, Acceleration
        REAL :: RIJ(3), RSQUARE, EPS=1.0, SIGMA=1.0,FORC
        REAL :: POT_EN, F(3), MASS=1.0, LAT,VIRIAL
      !
        E_POT = 0.0
        ACCN = 0.0
        VIRIAL = 0.0
      !
        DO I = 1, N-1
           DO J = I+1, N
              RIJ = ATOM_POS(I,:) - ATOM_POS(J,:)
              DO K = 1, 3
                 IF (ABS(RIJ(K)) > (L*LAT/2)) THEN
                    IF (RIJ(K) < 0) THEN  ! Sign flipped due to periodicity
                        RIJ(K) = L*LAT + RIJ(K)
                    ELSE
                        RIJ(K) = RIJ(K) - L*LAT
                    ENDIF
                 ENDIF
              ENDDO
      !
              RSQUARE = DOT_PRODUCT(RIJ, RIJ)
              POT_EN = 4.0*EPS*((SIGMA**12/RSQUARE**6)                  &
     &                         -(SIGMA**6/RSQUARE**3))
      !       V(Rij) = 4*epsilon*((sigma/Rij)^12-(sigma/Rij)^6)
              FORC = (48*EPS*(SIGMA**12/RSQUARE**7)                     &
     &           -24*EPS*(SIGMA**6/RSQUARE**4))   ! Force calculation
      !
              VIRIAL = VIRIAL - FORC * RSQUARE  
      !       F = -1/Rij(dV/dRij)
              E_POT(I) =  E_POT(I) + 0.5* POT_EN
              E_POT(J) =  E_POT(J) + 0.5* POT_EN
      !       Normalize Rij to find unit vectors
      !        RUIJ = RIJ/SQRT(RSQUARE)
      !       Update acceleration
              F = FORC * RIJ
              ACCN(I,:) = ACCN(I,:) + (F/MASS) ! Fij = -Fji
              ACCN(J,:) = ACCN(J,:) - (F/MASS)
           ENDDO
        ENDDO
        VIRIAL = - VIRIAL/3
      !
      END SUBROUTINE FORCE
