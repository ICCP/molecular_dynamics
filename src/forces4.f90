 
      subroutine calc_forces(N,pos,forces,pot,l,distances,presvirial)

      !implicit none
      integer, intent(in)  :: N
      real(8), intent(in)  :: pos(N,3),l
      real(8), intent(inout)  :: forces(N,3),pot, distances(N,N)
      real(8), intent(inout) :: presvirial(1)        
!f2py intent (in,out) :: forces, pot, distances, presvirial
      real(8), parameter :: rc =3.2
      real(8) :: dr_vec(3),partialforce(3),dr2,F,partialpot,cnt
      integer :: i,j
      
          
      presvirial = 0.0
      forces(:,:) = 0.0 
      pot = 0.0
      distances(:,:) = 0.0
      cnt = 0.0
      
      !calculating the radial distance between each pair of particles
      !Afterwards, do forces
      do j = 1, N
        do i = j+1, N
          dr_vec = pos(i,:) - pos(j,:)
          dr_vec = dr_vec - nint(dr_vec/l)*l
          dr2 = dot_product(dr_vec,dr_vec)
          distances(i,j) = sqrt(dr2)
          distances(j,i) = sqrt(dr2)
          if (dr2 < rc*rc) then
            dr2 = 1/dr2
            f = 2*dr2**7 - dr2**4
            partialforce(:) = 24*dr_vec*f
            forces(j,:) = forces(j,:) - partialforce(:)
            forces(i,:) = forces(i,:) + partialforce(:)
            presvirial = presvirial + dot_product(dr_vec,partialforce)
            partialpot = 4*(dr2**6 - dr2**3)
            pot = pot + partialpot
            cnt = cnt + 1
          end if
        enddo
      enddo
      presvirial = presvirial/cnt
       
      end subroutine calc_forces
      
    
