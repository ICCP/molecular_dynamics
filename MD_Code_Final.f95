Program MD_Code

IMPLICIT NONE
Integer :: i,count,count2,j
integer, parameter :: l=3
Integer, parameter :: N=4*(l**3)
Integer, parameter :: timeSteps=10000
real :: T0,dt,boxLength,time,kineticEnergy,potentialEnergy
real :: totalEnergy,B,T,P
Real, dimension(N,3) :: r0,r,v,F,F0
Real, dimension(3) :: totalF
integer, parameter :: out_unit1=10,out_unit2=20,out_unit3=30
integer, parameter :: out_unit4=40,out_unit5=50
!integer, parameter :: out_unit7=70,out_unit8=80,out_unit6=60
integer, parameter :: out_unit10=100,out_unit9=90
Real, dimension(timeSteps) :: MSD

T0=0.5
dt=0.0001
time=0

B=(2.0**(2.0/3.0))/2.0
boxLength=B*2*l

Call open_output_files()
Call seed_number_generator()

Call setup_lattice(l,r0,r,boxLength)
Call initialize_velocities(v,T0)

do i=1,N
    write (out_unit1,*) r(i,:)
    write (out_unit2,*) v(i,:)
end do

Call Calculate_Forces(r,F,boxLength,potentialEnergy,P)
Call Calculate_Kinetic_Energy(v,kineticEnergy)
write (out_unit3,*) time, kineticEnergy
write (out_unit4,*) time, potentialEnergy
	
do count=1,timeSteps
    time=time+dt
	Call Verlet(r,dt,v,F,boxLength)
    Call Shift_Velocity_Center(v)
    count2=count/100
	if (count2==1 .OR. count2==2 .OR. count2==3 .OR. count2==4) then
        Call Andersen_Thermostat(v,T0,T)
    End if		
	if (count2==5) then
        Call Andersen_Thermostat(v,T0,T)
    End if

!    Call Calculate_Potential_Energy(r,boxLength,potentialEnergy)
    Call Calculate_Kinetic_Energy(v,kineticEnergy)
    totalEnergy=potentialEnergy+kineticEnergy 
	Call Calculate_MSD(r,r0,MSD,count)
	
    write (out_unit3,*) time, kineticEnergy
    write (out_unit4,*) time, potentialEnergy
    write (out_unit5,*) time, totalEnergy
	write(out_unit9,*) time, MSD(count)
	write(out_unit10,*) time, P	
end do
	
Call close_output_files()

	
	T=0
	do i=1,N
	    do j=1,3
	        T=T+(v(i,j)**2)   
	    end do
	end do
    T=T/(3*(N-1))
	write(*,*) T

contains

Subroutine open_output_files()
    Implicit None

    open (unit=out_unit1,file="Initial_Lattice",action="write",status="replace")
    open (unit=out_unit2,file="Initial_Velocities",action="write",status="replace")
    open (unit=out_unit3,file="Kinetic_Energy",action="write",status="replace")
    open (unit=out_unit4,file="Potential_Energy",action="write",status="replace")
    open (unit=out_unit5,file="Total_Energy",action="write",status="replace")
!    open (unit=out_unit6,file="New_Positions",action="write",status="replace")
!    open (unit=out_unit7,file="New_Velocities",action="write",status="replace")
!    open (unit=out_unit8,file="Forces",action="write",status="replace")
    open (unit=out_unit9,file="MSD",action="write",status="replace")
	open (unit=out_unit10,file="Pressure",action="write",status="replace")
End Subroutine

Subroutine close_output_files()
    Implicit None

    close(out_unit1)
    close(out_unit2)
    close(out_unit3)
    close(out_unit4)
    close(out_unit5)
!    close(out_unit6)
!    close(out_unit7)
!    close(out_unit8)
	close(out_unit9)
	close(out_unit10)
End Subroutine	
	
Subroutine seed_number_generator()
    Implicit None
    integer :: i_seed
    integer, dimension(:), ALLOCATABLE :: Set_seed
    integer, dimension(1:8) :: dateSet_seed
    
	CALL RANDOM_SEED(size=i_seed)
    ALLOCATE(Set_seed(1:i_seed))
    CALL RANDOM_SEED(get=Set_seed)
    CALL DATE_AND_TIME(values=dateSet_seed)
    Set_seed(i_seed)=dateSet_seed(8)
    Set_seed(1)=dateSet_seed(8)*dateSet_seed(6)
    CALL RANDOM_SEED(put=Set_seed)
    DEALLOCATE(Set_seed)
End Subroutine

Subroutine setup_lattice(l,r0,r,boxLength)
    IMPLICIT NONE
    Real :: B,boxLength
    Integer :: l, i, count, j, k, m
    Real, dimension(:,:) :: r0,r
    Real, dimension(N,3) :: rUnit
 
    B=(2.0**(2.0/3.0))/2.0
	boxLength=B*2*l

    rUnit(:,:)=0.0
    rUnit(2,1)=B*1.0; rUnit(2,2)=B*1.0
    rUnit(3,1)=B*1.0; rUnit(3,3)=B*1.0
    rUnit(4,2)=B*1.0; rUnit(4,3)=B*1.0

count=1
do i=1,l
    do j=1,l
        do k=1,l
            do m=1,4
                r0(count,1)=rUnit(m,1)+Float(i-1)*B*2.0
                r0(count,2)=rUnit(m,2)+Float(j-1)*B*2.0
                r0(count,3)=rUnit(m,3)+Float(k-1)*B*2.0
				
                count=count+1
            end do
         end do
    end do
end do

do i=1,N
    r(i,:)=r0(i,:)
end do

End Subroutine

Subroutine get_random_number(u)
    IMPLICIT NONE
	
    real :: r,u

    CALL RANDOM_NUMBER(r)
    u = r
End Subroutine

Subroutine Shift_Velocity_Center(v)
    IMPLICIT NONE
    real, dimension(N,3) :: v
    real, dimension(3) :: sumV
    real, dimension(3) :: avgV
    integer :: i,j
	
    sumV(1)=0
    sumV(2)=0
    sumV(3)=0
    do i=1,N
        sumV(1)=sumV(1)+v(i,1)
        sumV(2)=sumV(2)+v(i,2)
        sumV(3)=sumV(3)+v(i,3)
    end do

        avgV(1)=sumV(1)/(N*1.0)
        avgV(2)=sumV(2)/(N*1.0)
        avgV(3)=sumV(3)/(N*1.0)
    do i=1,N
        do j=1,3
             v(i,j)=v(i,j)-avgV(j)
        end do	 
    end do
    
End Subroutine 

Subroutine initialize_velocities(v,T0)
    IMPLICIT NONE
    real, dimension(N,3) :: v
    real, dimension(3) :: sumV
    real :: u1,u2,u3,u4,T0
    real :: z1,z2,z3,z4,PI
    integer :: i
    PI=3.14159265359
    sumV(1)=0
    sumV(2)=0
    sumV(3)=0
    Call get_random_number(u1)
    !The first random number appears to always be 
    !0.15... so it is overwritten

    do i=1,N
        Call get_random_number(u1)
        Call get_random_number(u2)
        call get_random_number(u3)
        Call get_random_number(u4)
    
        z1=SQRT((-2)*log(u1))*cos((2*PI*u2))
        z2=SQRT((-2)*log(u1))*sin((2*PI*u2))
        z3=SQRT((-2)*log(u3))*cos((2*PI*u4))
        z4=SQRT((-2)*log(u3))*sin((2*PI*u4))
		
        v(i,1)=z1*(T0**0.5)
        v(i,2)=z2*(T0**0.5)
        v(i,3)=z3*(T0**0.5)
    

	end do

    Call Shift_Velocity_Center(v)
	
End Subroutine


Subroutine Calculate_Forces(r,F,boxLength,potentialEnergy,P)

    IMPLICIT NONE
	Real, dimension(N,3) :: r,F
    Real, dimension(3) :: dr,Ftemp
	Real :: boxLength,rsq,P
	Real :: potentialEnergy
	Integer :: i,j
	F(:,:)=0
	P=0.0
		
	potentialEnergy=0.0
	do i=1,N-1
		do j=i+1,N
		    if (i/=j) then
			   dr(:)=r(j,:)-r(i,:)
			   dr(:)=dr(:)-NINT(dr(:)/boxLength)*boxLength
			   rsq=dot_product(dr, dr)
               Ftemp(:)=(-48/(rsq**7) + 24/(rsq**4))*dr
			   potentialEnergy = potentialEnergy+(4.0*(1.0/rsq**6.0 - 1.0/rsq**3.0))
               F(i,:)=F(i,:)+Ftemp(:)
			   F(j,:)=F(j,:)-Ftemp(:)
			   P=P+dot_product(dr, Ftemp(:))                           
			end if
		end do
	end do
    P=P/3.0
	P=P+(N*T*1.0)
	P=P/((boxLength**3)*1.0)
End Subroutine

Subroutine Calculate_New_Positions(r,dt,v,F)
    IMPLICIT NONE
	Real, dimension(N,3) :: r,v,F,r0
	Real :: dt
	Integer :: i
	
	r0(:,:)=r(:,:)

	do i=1,N
        r(i,:)=r0(i,:)+(dt*v(i,:))+(0.5*(dt**2)*F(i,:))
    end do
End Subroutine

Subroutine Calculate_New_Velocities(v,dt,F,F0)
    IMPLICIT NONE
	Real, dimension(N,3) :: r,v,F,v0,F0
	Real :: dt
	Integer :: i
	
	v0(:,:)=v(:,:)
	do i=1,N
        v(i,:)=v0(i,:)+(dt*0.5*(F(i,:)+F0(i,:)))
    end do
End Subroutine

Subroutine Calculate_Kinetic_Energy(v,kineticEnergy)
    IMPLICIT NONE
	Real, dimension(N,3) :: v
	Real :: kineticEnergy
	Integer :: i
	
	kineticEnergy=0
	do i=1,N
	    kineticEnergy=kineticEnergy+(v(i,1)**2)+(v(i,2)**2)+(v(i,3)**2)
	end do
	kineticEnergy=kineticEnergy*0.5
End Subroutine

Subroutine Verlet(r,dt,v,F,boxLength)
    IMPLICIT NONE
	Real, dimension(N,3) :: r,v,F
	Real :: dt,boxLength
	
	v(:,:)=v(:,:)+(0.5*dt*F(:,:))
	r(:,:)=r(:,:)+(v(:,:)*dt)
	Call Calculate_Forces(r,F,boxLength,potentialEnergy,P)
	v(:,:)=v(:,:)+(0.5*F(:,:)*dt)
End Subroutine	
	
Subroutine Calculate_Temperature(v,T)
    IMPLICIT NONE
	Real, dimension(N,3) :: v
	Real :: T
	Integer :: i,j
	
	T=0
	do i=1,N
	    do j=1,3
	        T=T+(v(i,j)**2)   
	    end do
	end do
    T=T/(3*(N-1))	
End Subroutine
	
Subroutine Andersen_Thermostat(v,T0,T)
    IMPLICIT NONE
	Real, dimension(N,3) :: v
	Real :: T0,T
	
    Call Calculate_Temperature(v,T)
	v(:,:)=v(:,:)*sqrt((T0/T))

	
End Subroutine

Subroutine Calculate_MSD(r,r0,MSD,count)
    IMPLICIT NONE
	Real, dimension(N,3) :: r,r0
	Real, dimension(timeSteps) :: MSD
    Integer :: count
	MSD(count)=0.0
	
	do i=1,N
	    do j=1,3
		    MSD(count)=MSD(count)+((r(i,j)-r0(i,j))**2)
	    end do
	end do
	MSD(count)=MSD(count)/4.0

End Subroutine


end program MD_Code
