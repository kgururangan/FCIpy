module ci

        implicit none

        contains

       subroutine update_dpr(q,&
                             ndet,&
                             e_diagonal,&
                             energy)

                integer, intent(in) :: ndet
                double precision, intent(inout) :: q(ndet)
                double precision, intent(in) :: energy
                double precision, intent(in) :: e_diagonal(ndet)
                !
                integer :: idet
                double precision :: e0

                do idet=1,ndet
                   q(idet) = q(idet)/(energy - e_diagonal(idet) + 1.0d-012)
                end do

       end subroutine update_dpr

       subroutine calc_diagonal(det,N_int,ndet,&
                                noa,nob,&
                                e1int,e2int,mo_num,&
                                e_diag)

                integer, intent(in) :: N_int, ndet, mo_num, noa, nob
                integer*8, intent(in) :: det(N_int,2,ndet)
                double precision, intent(in) :: e1int(mo_num,mo_num)
                double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
                !
                double precision, intent(out) :: e_diag(ndet)
                !
                integer :: idet
                integer :: occ(noa,2)
                double precision :: e0

                do idet=1,ndet
                   call get_occupied(det(:,:,idet), N_int, noa, occ)
                   call slater0(occ,noa,nob,mo_num,e1int,e2int,e0)
                   e_diag(idet) = e0
                end do

       end subroutine calc_diagonal

       subroutine calc_sigma(det,N_int,ndet,&
                             e1int,e2int,mo_num,&
                             noa,nob,&
                             coef,&
                             sigma)

                integer, intent(in) :: N_int, ndet, mo_num, noa, nob
                integer*8, intent(in) :: det(N_int,2,ndet)
                double precision, intent(in) :: e1int(mo_num,mo_num)
                double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
                double precision, intent(in) :: coef(ndet)
                !
                double precision, intent(out) :: sigma(ndet)
                !
                integer :: idet, jdet
                integer :: occ(noa,2)
                double precision :: hmatel, phase
                integer :: exc(0:2,2,2)
                integer :: degree

                sigma = 0.0d0
                !
                ! Diagonal parts
                !
                do idet=1,ndet
                   call get_occupied(det(:,:,idet), N_int, noa, occ)
                   call slater0(occ,noa,nob,mo_num,e1int,e2int,hmatel)
                   sigma(idet) = hmatel*coef(idet)
                end do
                !
                ! Off-Diagonal parts
                !
                do idet=1,ndet
                   call get_occupied(det(:,:,idet), N_int, noa, occ)
                   do jdet=idet+1,ndet
                      hmatel = 0.0d0
                      call get_excitation(det(:,:,idet),det(:,:,jdet),exc,degree,phase,N_int)
                      
                      if (degree > 2) cycle
                      
                      if (degree==1) then
                         call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                      else if (degree==2) then
                         call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                      end if
                      
                      sigma(idet) = sigma(idet)+hmatel*coef(jdet)
                      sigma(jdet) = sigma(jdet)+hmatel*coef(idet)
                   end do
                end do
       end subroutine calc_sigma

       subroutine calc_sigma_nonhermitian(det,N_int,ndet,&
                                          e1int,e2int,mo_num,&
                                          noa,nob,&
                                          coef,&
                                          sigma)

                integer, intent(in) :: N_int, ndet, mo_num, noa, nob
                integer*8, intent(in) :: det(N_int,2,ndet)
                double precision, intent(in) :: e1int(mo_num,mo_num)
                double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
                double precision, intent(in) :: coef(ndet)
                !
                double precision, intent(out) :: sigma(ndet)
                !
                integer :: idet, jdet
                integer :: occ(noa,2)
                double precision :: hmatel, phase
                integer :: exc(0:2,2,2)
                integer :: degree

                sigma = 0.0d0
                do idet=1,ndet
                   call get_occupied(det(:,:,idet), N_int, noa, occ)
                   do jdet=1,ndet
                      hmatel = 0.0d0
                      call get_excitation(det(:,:,idet),det(:,:,jdet),exc,degree,phase,N_int)
                      
                      if (degree > 2) cycle

                      if (degree==0) then
                         call slater0(occ,noa,nob,mo_num,e1int,e2int,hmatel)
                      else if (degree==1) then
                         call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                      else if (degree==2) then
                         call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                      end if
                      
                      sigma(idet) = sigma(idet)+hmatel*coef(jdet)
                   end do
                end do
       end subroutine calc_sigma_nonhermitian
   
       subroutine calc_sigma_opt(det,N_int,ndet,num_alpha,num_beta,&
                                 e1int,e2int,mo_num,&
                                 noa,nob,&
                                 coef,&
                                 sigma,&
                                 e_diag)

                integer, intent(in) :: N_int, ndet, mo_num, noa, nob, num_alpha, num_beta
                double precision, intent(in) :: e1int(mo_num,mo_num)
                double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
                !
                integer*8, intent(inout) :: det(N_int,2,ndet)
                !f2py intent(in,out) :: det(N_int,2,ndet)
                double precision, intent(inout) :: coef(ndet)
                !f2py intent(in,out) :: coef(ndet)
                double precision, intent(inout) :: e_diag(ndet)
                !f2py intent(in,out) :: e_diag(ndet)

                double precision, intent(out) :: sigma(ndet)
                !
                integer :: a1, a2, b1, b2, a, b
                integer :: idet, jdet
                integer :: occ(noa,2)
                double precision :: hmatel, phase
                integer :: exc(0:2,2,2)
                integer :: degree
                integer, allocatable :: loc_arr(:)

                !
                ! diagonal loop
                !
                do idet=1,ndet
                   ! sigma(I) <- <I|H|I>*I
                   sigma(idet) = e_diag(idet)*coef(idet)
                end do
                !
                ! alpha-alpha loop
                !
                ! sort determinants into blocks according to beta string
                allocate(loc_arr(num_beta+1))
                call sort_determinants(sigma,coef,det,N_int,ndet,loc_arr,num_beta,2,e_diag)
                !
                do a=1,num_beta
                   do b1=loc_arr(a),loc_arr(a+1)-1
                      call get_occupied(det(:,:,b1), N_int, noa, occ)
                      do b2=b1+1,loc_arr(a+1)-1
                         if (exc_degree(det(:,1,b1),det(:,1,b2),N_int) <= 2) then
                            ! b1 connected to b2 by a or aa excitation
                            call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                            if (degree==1) then
                               call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               ! sigma(b1) <- <b1|H|b2>*c_b2
                               sigma(b1) = sigma(b1)+hmatel*coef(b2)
                               ! sigma(b2) <- <b2|H|b1>*c_b1
                               sigma(b2) = sigma(b2)+hmatel*coef(b1)
                            else if (degree==2) then
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               ! sigma(b1) <- <b1|H|b2>*c_b2
                               sigma(b1) = sigma(b1)+hmatel*coef(b2)
                               ! sigma(b2) <- <b2|H|b1>*c_b1
                               sigma(b2) = sigma(b2)+hmatel*coef(b1)
                            end if
                         end if
                      end do
                   end do
                end do
                !print*, 'Checking beta spin buckets'
                !print*, '---------------------------'
                !call check_spin_buckets(N_int, ndet, num_beta, loc_arr, det)
                deallocate(loc_arr)
                !
                ! beta-beta loop
                !
                ! sort determinants into blocks according to alpha string
                allocate(loc_arr(num_alpha+1))
                call sort_determinants(sigma,coef,det,N_int,ndet,loc_arr,num_alpha,1,e_diag)
                !
                do a=1,num_alpha
                   do b1=loc_arr(a),loc_arr(a+1)-1
                      call get_occupied(det(:,:,b1), N_int, noa, occ)
                      do b2=b1+1,loc_arr(a+1)-1
                         if (exc_degree(det(:,2,b1),det(:,2,b2),N_int) <= 2) then
                            ! b1 connected to b2 by b or bb excitation
                            call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                            if (degree==1) then
                               call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               ! sigma(b1) <- <b1|H|b2>*c_b2
                               sigma(b1) = sigma(b1)+hmatel*coef(b2)
                               ! sigma(b2) <- <b2|H|b1>*c_b1
                               sigma(b2) = sigma(b2)+hmatel*coef(b1)
                            else if (degree==2) then
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               ! sigma(b1) <- <b1|H|b2>*c_b2
                               sigma(b1) = sigma(b1)+hmatel*coef(b2)
                               ! sigma(b2) <- <b2|H|b1>*c_b1
                               sigma(b2) = sigma(b2)+hmatel*coef(b1)
                            end if
                         end if
                      end do
                   end do
                end do
                !print*, 'Checking alpha spin buckets'
                !print*, '---------------------------'
                !call check_spin_buckets(N_int, ndet, num_alpha, loc_arr, det)
                !
                ! alpha-beta loop
                !
                do a1=1,num_alpha
                   do a2=a1+1,num_alpha
                      ! get alpha excitation level between buckets
                      if (exc_degree(det(:,1,loc_arr(a1)),det(:,1,loc_arr(a2)),N_int) /= 1) cycle
                      do b1=loc_arr(a1),loc_arr(a1+1)-1
                         call get_occupied(det(:,:,b1), N_int, noa, occ)
                         do b2=loc_arr(a2),loc_arr(a2+1)-1
                            if (exc_degree(det(:,2,b1),det(:,2,b2),N_int) == 1) then
                               ! b1 connected to b2 by ab excitations
                               call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               ! sigma(b1) <- <b1|H|b2>*c_b2
                               sigma(b1) = sigma(b1)+hmatel*coef(b2)
                               ! sigma(b2) <- <b2|H|b1>*c_b1
                               sigma(b2) = sigma(b2)+hmatel*coef(b1)
                            end if
                         end do
                      end do
                   end do
                end do
                deallocate(loc_arr)
       end subroutine calc_sigma_opt
   
       subroutine build_hamiltonian(det,N_int,ndet,e1int,e2int,mo_num,noa,nob,hmat)

            integer, intent(in)    :: N_int
            integer, intent(in) :: ndet
            integer, intent(in) :: mo_num, noa, nob
            integer*8, intent(in)  :: det(N_int,2,ndet)
            double precision, intent(in) :: e1int(mo_num,mo_num)
            double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
            !
            double precision, intent(out) :: hmat(ndet,ndet)
            !
            integer          :: occ(noa,2)
            double precision :: zero, one, two
            double precision :: phase
            integer    :: exc(0:2,2,2)
            integer    :: degree
            integer :: i, j

            hmat=0.0
            do i=1,ndet
              call get_occupied(det(:,:,i), N_int, noa, occ)
              do j=i,ndet
                zero = 0.0
                one = 0.0
                two = 0.0
                call get_excitation(det(:,:,i),det(:,:,j),exc,degree,phase,N_int)
                if (degree==0) then
                  call slater0(occ,noa,nob,mo_num,e1int,e2int,zero)
                else if (degree==1) then
                  call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,one)
                else if (degree==2) then
                  call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,two)
                end if
                hmat(i,j)=zero+one+two
                hmat(j,i)=hmat(i,j)
              end do
            end do

       end subroutine build_hamiltonian
   
       subroutine build_hamiltonian_opt(det,N_int,ndet,num_alpha,num_beta,&
                                        e1int,e2int,mo_num,&
                                        noa,nob,&
                                        H)

                integer, intent(in) :: N_int, ndet, mo_num, noa, nob, num_alpha, num_beta
                double precision, intent(in) :: e1int(mo_num,mo_num)
                double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
                !
                integer*8, intent(inout) :: det(N_int,2,ndet)
                !f2py intent(in,out) :: det(N_int,2,ndet)
                double precision, intent(out) :: H(ndet,ndet)
                !
                integer :: a1, a2, b1, b2, a, b
                integer :: idet, jdet
                integer :: occ(noa,2)
                double precision :: hmatel, phase
                integer :: exc(0:2,2,2)
                integer :: degree
                integer, allocatable :: loc_arr(:)
                !integer*8 :: det1, det2

                !
                ! diagonal loop
                !
                do idet=1,ndet
                   hmatel = 0.0d0
                   call get_occupied(det(:,:,idet), N_int, noa, occ)
                   call slater0(occ,noa,nob,mo_num,e1int,e2int,hmatel)
                   H(idet,idet) = hmatel
                end do
                !
                ! alpha-alpha loop
                !
                ! sort determinants into blocks according to beta string
                allocate(loc_arr(num_beta+1))
                call sort_determinants_simple(det,H,N_int,ndet,loc_arr,num_beta,2)
                !
                do a=1,num_beta
                   do b1=loc_arr(a),loc_arr(a+1)-1
                      call get_occupied(det(:,:,b1), N_int, noa, occ)
                      do b2=b1+1,loc_arr(a+1)-1
                         if (exc_degree(det(:,1,b1),det(:,1,b2),N_int) <= 2) then
                            hmatel = 0.0d0
                            ! b1 connected to b2 by a or aa excitation
                            call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                            if (degree==1) then
                               call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                            else if (degree==2) then
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                            end if
                            H(b1,b2) = hmatel
                            H(b2,b1) = H(b1,b2)
                         end if
                      end do
                   end do
                end do
                !print*, 'Checking beta spin buckets'
                !print*, '---------------------------'
                !call check_spin_buckets(N_int, ndet, num_beta, loc_arr, det)
                deallocate(loc_arr)
                !
                ! beta-beta loop
                !
                ! sort determinants into blocks according to alpha string
                allocate(loc_arr(num_alpha+1))
                call sort_determinants_simple(det,H,N_int,ndet,loc_arr,num_alpha,1)
                !
                do a=1,num_alpha
                   do b1=loc_arr(a),loc_arr(a+1)-1
                      call get_occupied(det(:,:,b1), N_int, noa, occ)
                      do b2=b1+1,loc_arr(a+1)-1
                         if (exc_degree(det(:,2,b1),det(:,2,b2),N_int) <= 2) then
                            hmatel = 0.0d0
                            ! b1 connected to b2 by b or bb excitation
                            call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                            if (degree==1) then
                               call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                            else if (degree==2) then
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                            end if
                            H(b1,b2) = hmatel
                            H(b2,b1) = H(b1,b2)
                         end if
                      end do
                   end do
                end do
                !print*, 'Checking alpha spin buckets'
                !print*, '---------------------------'
                !call check_spin_buckets(N_int, ndet, num_alpha, loc_arr, det)
                !
                ! alpha-beta loop
                !
                do a1=1,num_alpha
                   do a2=a1+1,num_alpha
                      ! get alpha excitation level between buckets
                      if (exc_degree(det(:,1,loc_arr(a1)),det(:,1,loc_arr(a2)),N_int) /= 1) cycle
                      do b1=loc_arr(a1),loc_arr(a1+1)-1
                         call get_occupied(det(:,:,b1), N_int, noa, occ)
                         do b2=loc_arr(a2),loc_arr(a2+1)-1
                            if (exc_degree(det(:,2,b1),det(:,2,b2),N_int) == 1) then
                               hmatel = 0.0d0
                               ! b1 connected to b2 by ab excitations
                               call get_excitation(det(:,:,b1),det(:,:,b2),exc,degree,phase,N_int)
                               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                               H(b1,b2) = hmatel
                               H(b2,b1) = H(b1,b2)
                            end if
                         end do
                      end do
                   end do
                end do
                deallocate(loc_arr)

       end subroutine build_hamiltonian_opt
   
       subroutine build_hamiltonian_nonhermitian(det,N_int,ndet,e1int,e2int,mo_num,noa,nob,hmat)

            integer, intent(in)    :: N_int
            integer, intent(in) :: ndet
            integer, intent(in) :: mo_num, noa, nob
            integer*8, intent(in)  :: det(N_int,2,ndet)
            double precision, intent(in) :: e1int(mo_num,mo_num)
            double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
            !
            double precision, intent(out) :: hmat(ndet,ndet)
            !
            integer          :: occ(noa,2)
            double precision :: zero, one, two
            double precision :: phase
            integer    :: exc(0:2,2,2)
            integer    :: degree
            integer :: i, j

            hmat=0.0
            do i=1,ndet
              call get_occupied(det(:,:,i), N_int, noa, occ)
              do j=1,ndet
                zero = 0.0
                one = 0.0
                two = 0.0
                call get_excitation(det(:,:,i),det(:,:,j),exc,degree,phase,N_int)
                if (degree==0) then
                  call slater0(occ,noa,nob,mo_num,e1int,e2int,zero)
                else if (degree==1) then
                  call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,one)
                else if (degree==2) then
                  call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,two)
                end if
                hmat(i,j)=zero+one+two
              end do
            end do

       end subroutine build_hamiltonian_nonhermitian

       subroutine i_H_j(I,J,N_int,e1int,e2int,mo_num,noa,nob,hmatel)
            integer, intent(in)    :: N_int
            integer, intent(in) :: mo_num, noa, nob
            integer*8, intent(in)  :: I(N_int,2), J(N_int,2)
            double precision, intent(in) :: e1int(mo_num,mo_num)
            double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)
            !
            double precision, intent(out) :: hmatel
            !
            integer          :: occ(noa,2)
            double precision :: zero, one, two
            double precision :: phase
            integer    :: exc(0:2,2,2)
            integer    :: degree
           
            hmatel = 0.0

            call get_occupied(I, N_int, noa, occ)
            call get_excitation(I,J,exc,degree,phase,N_int)
            if (degree==0) then
               call slater0(occ,noa,nob,mo_num,e1int,e2int,hmatel)
            else if (degree==1) then
               call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
            else if (degree==2) then
               call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
            end if
       end subroutine i_H_j

       subroutine get_occupied(det,Nint,noa,occ)
          implicit none
          integer, intent(in) :: Nint
          integer*8, intent(in) :: det(Nint,2)
          integer, intent(in)   :: noa
          integer, intent(out)  :: occ(noa,2)

          integer :: ispin, ishift, i, j
          integer*8 :: buffer
          integer :: n
         
          do ispin=1,2
            n = 1
            ishift = 0
            do i=1,Nint
              buffer = det(i,ispin)
              do while (buffer /= 0_8)
                j = trailz(buffer) + ishift
                buffer = iand(buffer,buffer-1_8)
                occ(n,ispin) = j+1
                n = n+1
              end do
              ishift = ishift+64
            end do
          end do
        end 

        subroutine slater0(occ,noa,nob,mo_num,e1int,e2int,val)
          implicit none
          integer, intent(in)   :: noa, nob, mo_num
          integer, intent(in)  :: occ(noa,2)
          double precision, intent(in) :: e1int(mo_num,mo_num)
          double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)

          integer :: i, j, k, l
          double precision :: val

          val = 0.0
          ! onebody a
          do i=1,noa
            k=occ(i,1)
            val = val+e1int(k,k)
          end do
          ! onebody b
          do i=1,nob
            k=occ(i,2)
            val = val+e1int(k,k)
          end do
          ! twobody aa
          do i=1,noa
            k=occ(i,1)
            do j=i+1,noa
              l=occ(j,1)
              val = val+(e2int(k,l,k,l)-e2int(k,l,l,k))
            end do
          end do
          ! twobody ab
          do i=1,noa
            k=occ(i,1)
            do j=1,nob
              l=occ(j,2)
              val = val+e2int(k,l,k,l)
            end do
          end do
          ! twobody bb
          do i=1,nob
            k=occ(i,2)
            do j=i+1,nob
              l=occ(j,2)
              val = val+(e2int(k,l,k,l)-e2int(k,l,l,k))
            end do
          end do
        end

        subroutine slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,val)
          implicit none
          integer, intent(in)   :: noa, nob, mo_num
          integer, intent(in)  :: occ(noa,2)
          integer, intent(out) :: exc(0:2,2,2)
          double precision, intent(in) :: phase
          double precision, intent(in) :: e1int(mo_num,mo_num)
          double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)

          integer :: i, j, k, l
          double precision :: val
          logical :: alpha

          val = 0.0
          if (exc(0,1,1)==1) then ! single excitation is a->a
             alpha=.true.
             i=exc(1,1,1)
             j=exc(1,2,1)
          else ! single excitation is b->b
             alpha=.false.
             i=exc(1,1,2)
             j=exc(1,2,2)
          end if
          val = e1int(i,j)

          if (alpha) then
             ! twobody aa
             do k=1,noa
                l=occ(k,1)
                val = val + e2int(i,l,j,l) - e2int(i,l,l,j)
             end do
             ! twobody ab
             do k=1,nob
                l=occ(k,2)
                val = val + e2int(i,l,j,l)
             end do
          else
             ! twobody bb
             do k=1,nob
                l=occ(k,2)
                val = val + e2int(i,l,j,l) - e2int(i,l,l,j)
             end do
             ! twobody ab
             do k=1,noa
                l=occ(k,1)
                val = val + e2int(l,i,l,j)
             end do
          end if
          val = val * phase
        end 

        subroutine slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,val)
          implicit none
          integer, intent(in)   :: noa, nob, mo_num
          integer, intent(in)  :: occ(noa,2)
          integer, intent(out) :: exc(0:2,2,2)
          double precision, intent(in) :: phase
          double precision, intent(in) :: e1int(mo_num,mo_num)
          double precision, intent(in) :: e2int(mo_num,mo_num,mo_num,mo_num)

          integer :: i, j, k, l
          double precision :: val
          logical :: aa, ab, bb

          aa = .false.
          bb = .false.
          ab = .false.

          val = 0.0
          if (exc(0,1,1)==2) then ! double excitation is aa->aa
             aa = .true.
             i=exc(1,1,1)
             j=exc(1,2,1)
             k=exc(2,1,1)
             l=exc(2,2,1)
          else if (exc(0,1,2)==2) then ! double excitation is bb->bb
             bb = .true.
             i=exc(1,1,2)
             j=exc(1,2,2)
             k=exc(2,1,2)
             l=exc(2,2,2)
          else ! double excitation is ab->ab
             ab = .true.
             i=exc(1,1,1)
             j=exc(1,2,1)
             k=exc(1,1,2)
             l=exc(1,2,2)
          end if
          if (aa .or. bb) then
            val = (e2int(i,k,j,l) - e2int(i,k,l,j)) * phase
          else
            val = e2int(i,k,j,l) * phase
          end if
        end 

        subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
         !   Slater-Condon Rules
         !   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
         !                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
         !
         !   This program is free software; you can redistribute it and/or modify
         !   it under the terms of the GNU General Public License as published by
         !   the Free Software Foundation; either version 2 of the License, or
         !   (at your option) any later version.
         !
         !   This program is distributed in the hope that it will be useful,
         !   but WITHOUT ANY WARRANTY; without even the implied warranty of
         !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         !   GNU General Public License for more details.
         !
         !   You should have received a copy of the GNU General Public License along
         !   with this program; if not, write to the Free Software Foundation, Inc.,
         !   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
         integer, intent(in)  :: Nint
         integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
         integer, intent(out) :: exc(0:2,2,2)
         integer, intent(out) :: degree
         double precision, intent(out) :: phase

         degree = n_excitations(det1,det2,Nint)

         select case (degree)
         
           case (3:)
             degree = -1
             return

           case (2)
             call get_double_excitation(det1,det2,exc,phase,Nint)
             return

           case (1)
             call get_single_excitation(det1,det2,exc,phase,Nint)
             return

           case(0)
             return

           end select
        end

        subroutine get_excitation3(det1,det2,exc,degree,phase,Nint)
         integer, intent(in)  :: Nint
         integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
         integer, intent(out) :: exc(0:3,2,2)
         integer, intent(out) :: degree
         double precision, intent(out) :: phase

         degree = n_excitations(det1,det2,Nint)

         select case (degree)
         
           case (4:)
             degree = -1
             return

           case (3)
             call get_triple_excitation(det1,det2,exc,phase,Nint)
             return

           case (2)
             call get_double_excitation(det1,det2,exc(:2,:,:),phase,Nint)
             return

           case (1)
             call get_single_excitation(det1,det2,exc(:2,:,:),phase,Nint)
             return

           case(0)
             return

           end select
        end


        subroutine get_single_excitation(det1,det2,exc,phase,Nint)
         !   Slater-Condon Rules
         !   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
         !                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
         !
         !   This program is free software; you can redistribute it and/or modify
         !   it under the terms of the GNU General Public License as published by
         !   the Free Software Foundation; either version 2 of the License, or
         !   (at your option) any later version.
         !
         !   This program is distributed in the hope that it will be useful,
         !   but WITHOUT ANY WARRANTY; without even the implied warranty of
         !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         !   GNU General Public License for more details.
         !
         !   You should have received a copy of the GNU General Public License along
         !   with this program; if not, write to the Free Software Foundation, Inc.,
         !   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
         integer, intent(in)  :: Nint
         integer*8, intent(in)  :: det1(Nint,2)
         integer*8, intent(in)  :: det2(Nint,2)
         integer, intent(out) :: exc(0:2,2,2)
         double precision, intent(out) :: phase
         integer :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
         integer*8 :: hole, particle, tmp
         double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

         exc(0,1,1) = 0
         exc(0,2,1) = 0
         exc(0,1,2) = 0
         exc(0,2,2) = 0
         do ispin = 1,2
           ishift = -63
           do l=1,Nint
             ishift = ishift + 64
             if (det1(l,ispin) == det2(l,ispin)) cycle
             tmp = xor( det1(l,ispin), det2(l,ispin) )
             particle = iand(tmp, det2(l,ispin))
             hole     = iand(tmp, det1(l,ispin))
             if (particle /= 0_8) then
               tz = trailz(particle)
               exc(0,2,ispin) = 1
               exc(1,2,ispin) = tz+ishift
             end if
             if (hole /= 0_8) then
               tz = trailz(hole)
               exc(0,1,ispin) = 1
               exc(1,1,ispin) = tz+ishift
             end if

             if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
               low  = min(exc(1,1,ispin),exc(1,2,ispin))
               high = max(exc(1,1,ispin),exc(1,2,ispin))
               j = ishft(low-1,-6)+1
               n = iand(low-1,63)
               k = ishft(high-1,-6)+1
               m = iand(high-1,63)
               if (j==k) then
                 nperm = popcnt(iand(det1(j,ispin), &
                    iand( not(ishft(1_8,n+1))+1_8 ,ishft(1_8,m)-1_8)))
               else
                 nperm = popcnt(iand(det1(k,ispin), ishft(1_8,m)-1_8)) + &
                         popcnt(iand(det1(j,ispin), not(ishft(1_8,n+1))+1_8))
                 do i=j+1,k-1
                   nperm = nperm + popcnt(det1(i,ispin))
                 end do
               end if
               phase = phase_dble(iand(nperm,1))
               return
             end if
           end do
         end do
        end

        subroutine get_double_excitation(det1,det2,exc,phase,Nint)
         !   Slater-Condon Rules
         !   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
         !                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
         !
         !   This program is free software; you can redistribute it and/or modify
         !   it under the terms of the GNU General Public License as published by
         !   the Free Software Foundation; either version 2 of the License, or
         !   (at your option) any later version.
         !
         !   This program is distributed in the hope that it will be useful,
         !   but WITHOUT ANY WARRANTY; without even the implied warranty of
         !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         !   GNU General Public License for more details.
         !
         !   You should have received a copy of the GNU General Public License along
         !   with this program; if not, write to the Free Software Foundation, Inc.,
         !   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
         integer, intent(in)  :: Nint
         integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
         integer, intent(out) :: exc(0:2,2,2)
         double precision, intent(out) :: phase
         integer :: l, ispin, idx_hole, idx_particle, ishift
         integer :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
         integer*8 :: hole, particle, tmp
         double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
         exc(0,1,1) = 0
         exc(0,2,1) = 0
         exc(0,1,2) = 0
         exc(0,2,2) = 0
         nexc=0
         nperm=0
         do ispin = 1,2
          idx_particle = 0
          idx_hole = 0
          ishift = -63
          do l=1,Nint
           ishift = ishift + 64
           if (det1(l,ispin) == det2(l,ispin))  then
             cycle
           end if
           tmp = xor( det1(l,ispin), det2(l,ispin) )
           particle = iand(tmp, det2(l,ispin))
           hole     = iand(tmp, det1(l,ispin))
           do while (particle /= 0_8)
             tz = trailz(particle)
             nexc = nexc+1
             idx_particle = idx_particle + 1
             exc(0,2,ispin) = exc(0,2,ispin) + 1
             exc(idx_particle,2,ispin) = tz+ishift
             particle = iand(particle,particle-1_8)
           end do
           do while (hole /= 0_8)
             tz = trailz(hole)
             nexc = nexc+1
             idx_hole = idx_hole + 1
             exc(0,1,ispin) = exc(0,1,ispin) + 1
             exc(idx_hole,1,ispin) = tz+ishift
             hole = iand(hole,hole-1_8)
           end do
           if (nexc == 4) exit
          end do

          do i=1,exc(0,1,ispin)
            low  = min(exc(i,1,ispin),exc(i,2,ispin))
            high = max(exc(i,1,ispin),exc(i,2,ispin))
            j = ishft(low-1,-6)+1
            n = iand(low-1,63)
            k = ishft(high-1,-6)+1
            m = iand(high-1,63)
            if (j==k) then
              nperm = nperm + popcnt(iand(det1(j,ispin),  &
                 iand( not(ishft(1_8,n+1))+1 ,ishft(1_8,m)-1)))
            else
              nperm = nperm + popcnt(iand(det1(k,ispin),  &
                                     ishft(1_8,m)-1)) &
                            + popcnt(iand(det1(j,ispin),  &
                                     not(ishft(1_8,n+1))+1))
              do l=j+1,k-1
                nperm = nperm + popcnt(det1(l,ispin))
              end do
            end if
          end do
          if (exc(0,1,ispin) == 2) then
            a = min(exc(1,1,ispin), exc(1,2,ispin))
            b = max(exc(1,1,ispin), exc(1,2,ispin))
            c = min(exc(2,1,ispin), exc(2,2,ispin))
            d = max(exc(2,1,ispin), exc(2,2,ispin))
            if (c>a .and. c<b .and. d>b) nperm = nperm + 1
            exit
           end if
         end do
         phase = phase_dble(iand(nperm,1))
        end

        subroutine get_triple_excitation(det1, det2, exc, phase, Nint)
        !
        !  Slater-Condon Rules: triple excitation finder
        !
              integer, intent(in)          :: Nint
              integer*8, intent(in)        :: det1(Nint,2), det2(Nint,2)
              integer, intent(out)         :: exc(0:3,2,2)
              double precision, intent(out):: phase

              integer*8 :: tmp, hole, particle, mask, tmp_lower
              integer    :: ispin, l, ishift, tz
              integer    :: nexc, idx_hole, idx_particle
              integer    :: i_h, i_p, segment, pos, seg, idx
              integer    :: nperm, parity
              integer*8 :: tempDet(Nint)
              double precision, parameter :: phase_table(0:1) = (/ 1.d0, -1.d0 /)

        !-- Initialize counts
              exc(0,1,1) = 0; exc(0,2,1) = 0
              exc(0,1,2) = 0; exc(0,2,2) = 0
              nexc = 0

        !-- Detect holes and particles until 3 excitations (6 bit changes)
              do ispin = 1, 2
                idx_hole     = 0
                idx_particle = 0
                ishift       = -63
                do l = 1, Nint
                  ishift = ishift + 64
                  if (det1(l,ispin) == det2(l,ispin)) cycle
                  tmp = ieor(det1(l,ispin), det2(l,ispin))

                  ! Particles: bits present in det2 but not det1
                  particle = iand(tmp, det2(l,ispin))
                  do while (particle /= 0_8)
                    tz = trailz(particle)
                    idx_particle = idx_particle + 1
                    exc(0,2,ispin) = exc(0,2,ispin) + 1
                    exc(idx_particle,2,ispin) = tz + ishift
                    particle = iand(particle, particle - 1_8)
                    nexc = nexc + 1
                  end do

                  ! Holes: bits present in det1 but not det2
                  hole = iand(tmp, det1(l,ispin))
                  do while (hole /= 0_8)
                    tz = trailz(hole)
                    idx_hole = idx_hole + 1
                    exc(0,1,ispin) = exc(0,1,ispin) + 1
                    exc(idx_hole,1,ispin) = tz + ishift
                    hole = iand(hole, hole - 1_8)
                    nexc = nexc + 1
                  end do

                  if (nexc == 6) exit
                end do
              end do

        !-- Compute Fermionic phase by sequential annihilation then creation
              nperm = 0
              do ispin = 1, 2
                ! Copy original bits
                tempDet = det1(:,ispin)

                ! Annihilate holes
                do i_h = 1, exc(0,1,ispin)
                  idx = exc(i_h,1,ispin)
                  segment = (idx - 1) / 64 + 1
                  pos = mod(idx - 1, 64)
                  mask = shiftl(1_8, pos)
                  tmp_lower = iand(tempDet(segment), mask - 1_8)
                  parity = popcnt(tmp_lower)
                  do seg = 1, segment - 1
                    parity = parity + popcnt(tempDet(seg))
                  end do
                  nperm = nperm + parity
                  tempDet(segment) = iand(tempDet(segment), not(mask))
                end do

                ! Create particles
                do i_p = 1, exc(0,2,ispin)
                  idx = exc(i_p,2,ispin)
                  segment = (idx - 1) / 64 + 1
                  pos = mod(idx - 1, 64)
                  mask = shiftl(1_8, pos)
                  tmp_lower = iand(tempDet(segment), mask - 1_8)
                  parity = popcnt(tmp_lower)
                  do seg = 1, segment - 1
                    parity = parity + popcnt(tempDet(seg))
                  end do
                  nperm = nperm + parity
                  tempDet(segment) = ior(tempDet(segment), mask)
                end do
              end do

              !-- Final phase
              !-- Flip the sign - f**k you, o4-mini-high
              phase = -1.0 * phase_table(mod(nperm,2))

        end subroutine get_triple_excitation


        integer function exc_degree(det1,det2,Nint)
          !   Slater-Condon Rules
          !   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
          !                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
          !
          !   This program is free software; you can redistribute it and/or modify
          !   it under the terms of the GNU General Public License as published by
          !   the Free Software Foundation; either version 2 of the License, or
          !   (at your option) any later version.
          !
          !   This program is distributed in the hope that it will be useful,
          !   but WITHOUT ANY WARRANTY; without even the implied warranty of
          !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
          !   GNU General Public License for more details.
          !
          !   You should have received a copy of the GNU General Public License along
          !   with this program; if not, write to the Free Software Foundation, Inc.,
          !   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
          integer*8, intent(in)  :: det1(Nint), det2(Nint)
          integer  , intent(in)  :: Nint

          integer                :: l
         
          exc_degree = &
              popcnt(xor( det1(1), det2(1)) )
         
          do l=2,Nint
            exc_degree = exc_degree + &
              popcnt(xor( det1(l), det2(l)) )
          end do
          exc_degree = ishft(exc_degree,-1)

        end


        integer function n_excitations(det1,det2,Nint)
          !   Slater-Condon Rules
          !   Copyright (C) 2013 Anthony Scemama <scemama@irsamc.ups-tlse.fr>
          !                      Emmanuel Giner <emmanuel_giner_jr@hotmail.fr>
          !
          !   This program is free software; you can redistribute it and/or modify
          !   it under the terms of the GNU General Public License as published by
          !   the Free Software Foundation; either version 2 of the License, or
          !   (at your option) any later version.
          !
          !   This program is distributed in the hope that it will be useful,
          !   but WITHOUT ANY WARRANTY; without even the implied warranty of
          !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
          !   GNU General Public License for more details.
          !
          !   You should have received a copy of the GNU General Public License along
          !   with this program; if not, write to the Free Software Foundation, Inc.,
          !   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
          integer*8, intent(in)  :: det1(Nint,2), det2(Nint,2)
          integer  , intent(in)  :: Nint

          integer                :: l
         
          n_excitations = &
              popcnt(xor( det1(1,1), det2(1,1)) ) + &
              popcnt(xor( det1(1,2), det2(1,2)) )
         
          do l=2,Nint
            n_excitations = n_excitations + &
              popcnt(xor( det1(l,1), det2(l,1)) ) + &
              popcnt(xor( det1(l,2), det2(l,2)) )
          end do
          n_excitations = ishft(n_excitations,-1)

        end

        subroutine compute_rdm1s(det,Ndet,coef,mo_num,noa,nob,&
                                 Nint,rdm1a,rdm1b)
         integer*8, intent(in)         :: det(Nint,2,Ndet)
         integer, intent(in)           :: Ndet, Nint, mo_num, noa, nob
         double precision, intent(in)  :: coef(Ndet)
         double precision, intent(out) :: rdm1a(mo_num,mo_num)
         double precision, intent(out) :: rdm1b(mo_num,mo_num)

         integer :: i,j,k,l,ispin,ishift
         integer*8 :: buffer, det_I(Nint,2), det_J(Nint,2)
         integer :: deg
         integer :: exc(0:2,2,2) ! [(0:n,1,2:mo_idx), hole/part, ispin]
         double precision :: phase, c
         integer :: occ(noa,2)

         rdm1a = 0.d0
         rdm1b = 0.d0
         do k=1,Ndet
          det_I = det(:,:,k)
          c = coef(k)*coef(k)
          call get_occupied(det_I,Nint,noa,occ)
          !
          ! diagonal: <I|p+ q|I>
          !
          ! alpha --->
          do j=1,noa
             rdm1a(occ(j,1),occ(j,1)) = rdm1a(occ(j,1),occ(j,1)) + c
          end do
          ! beta ---> 
          do j=1,nob
             rdm1b(occ(j,2),occ(j,2)) = rdm1b(occ(j,2),occ(j,2)) + c
          end do
          !
          ! off-diagonal: <J|p+ q|I>
          !
          do l=1,k-1
           det_J = det(:,:,l)
           call get_excitation(det_I,det_J,exc,deg,phase,Nint)
           if (deg /= 1) cycle
           c = phase*coef(k)*coef(l)
           if (exc(0,1,1) == 1) then ! a -> a
             i = exc(1,1,1)
             j = exc(1,2,1)
             rdm1a(j,i) = rdm1a(j,i) + c
             rdm1a(i,j) = rdm1a(i,j) + c
           else ! b -> b
             i = exc(1,1,2)
             j = exc(1,2,2)
             rdm1b(j,i) = rdm1b(j,i) + c
             rdm1b(i,j) = rdm1b(i,j) + c
           end if
          end do
         end do
        end

        subroutine compute_rdm2s(det,Ndet,coef,mo_num,noa,nob,&
                                 Nint,rdm2a,rdm2b,rdm2c)
         integer*8, intent(in)         :: det(Nint,2,Ndet)
         integer, intent(in)           :: Ndet, Nint, mo_num, noa, nob
         double precision, intent(in)  :: coef(Ndet)
         double precision, intent(out) :: rdm2a(mo_num,mo_num,mo_num,mo_num)
         double precision, intent(out) :: rdm2b(mo_num,mo_num,mo_num,mo_num)
         double precision, intent(out) :: rdm2c(mo_num,mo_num,mo_num,mo_num)

         integer :: i,j,k,l,m,n,ispin,ishift,ii,jj,kk,ll,mm,nn
         integer :: degree
         integer :: exc(0:2,2,2) ! [(0:n,1,2:mo_idx), hole/part, ispin]
         double precision :: phase, c, cp
         integer :: occ(noa,2)

         rdm2a = 0.d0
         rdm2b = 0.d0
         rdm2c = 0.d0
         do n=1,Ndet
            call get_occupied(det(:,:,n), Nint, noa, occ)
            c = coef(n)*coef(n)
            !
            ! Diagonal
            !
            ! aa --->
            do ii=1,noa
               do jj=ii+1,noa
                  i = occ(ii,1); j = occ(jj,1);
                  rdm2a(i,j,i,j) = rdm2a(i,j,i,j) + c
                  rdm2a(i,j,j,i) = rdm2a(i,j,j,i) - c
                  rdm2a(j,i,i,j) = rdm2a(j,i,i,j) - c
                  rdm2a(j,i,j,i) = rdm2a(j,i,j,i) + c
               end do
            end do
            ! ab ---> 
            do ii=1,noa
               do jj=1,nob
                  i = occ(ii,1); j = occ(jj,2);
                  rdm2b(i,j,i,j) = rdm2b(i,j,i,j) + c
               end do
            end do
            ! bb --->
            do ii=1,nob
               do jj=ii+1,nob
                  i = occ(ii,2); j = occ(jj,2);
                  rdm2c(i,j,i,j) = rdm2c(i,j,i,j) + c
                  rdm2c(i,j,j,i) = rdm2c(i,j,j,i) - c
                  rdm2c(j,i,i,j) = rdm2c(j,i,i,j) - c
                  rdm2c(j,i,j,i) = rdm2c(j,i,j,i) + c
               end do
            end do
            !
            ! Off-diagonal
            !
            do m=1,Ndet
               call get_excitation(det(:,:,n),det(:,:,m),exc,degree,phase,Nint)
               if (degree >= 3 .or. degree==0) cycle
               c = coef(n)*coef(m)
               cp = phase*c
               if (degree==1) then ! single excitation
                  if (exc(0,1,1)==1) then ! single excitation is a->a
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     ! aa --->
                     do kk=1,noa
                        k = occ(kk,1)
                        rdm2a(i,k,j,k) = rdm2a(i,k,j,k) + cp
                        rdm2a(k,i,j,k) = rdm2a(k,i,j,k) - cp
                        rdm2a(i,k,k,j) = rdm2a(i,k,k,j) - cp
                        rdm2a(k,i,k,j) = rdm2a(k,i,k,j) + cp
                     end do
                     ! ab ---> 
                     do kk=1,nob
                        k = occ(kk,2)
                        rdm2b(i,k,j,k) = rdm2b(i,k,j,k) + cp
                     end do
                  else ! single excitation is b->b
                     i=exc(1,1,2)
                     j=exc(1,2,2)
                     ! ab ---> 
                     do kk=1,noa
                        k = occ(kk,1)
                        rdm2b(k,i,k,j) = rdm2b(k,i,k,j) + cp
                     end do
                     ! bb --->
                     do kk=1,nob
                        k = occ(kk,2)
                        rdm2c(i,k,j,k) = rdm2c(i,k,j,k) + cp
                        rdm2c(k,i,j,k) = rdm2c(k,i,j,k) - cp
                        rdm2c(i,k,k,j) = rdm2c(i,k,k,j) - cp
                        rdm2c(k,i,k,j) = rdm2c(k,i,k,j) + cp
                     end do
                  end if
               else if (degree==2) then
                  if (exc(0,1,1)==2) then ! double excitation is aa->aa
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(2,1,1)
                     l=exc(2,2,1)
                     rdm2a(i,k,j,l) = rdm2a(i,k,j,l) + cp
                     rdm2a(k,i,j,l) = rdm2a(k,i,j,l) - cp
                     rdm2a(i,k,l,j) = rdm2a(i,k,l,j) - cp
                     rdm2a(k,i,l,j) = rdm2a(k,i,l,j) + cp
                  else if (exc(0,1,2)==2) then ! double excitation is bb->bb
                     i=exc(1,1,2)
                     j=exc(1,2,2)
                     k=exc(2,1,2)
                     l=exc(2,2,2)
                     rdm2c(i,k,j,l) = rdm2c(i,k,j,l) + cp
                     rdm2c(k,i,j,l) = rdm2c(k,i,j,l) - cp
                     rdm2c(i,k,l,j) = rdm2c(i,k,l,j) - cp
                     rdm2c(k,i,l,j) = rdm2c(k,i,l,j) + cp
                  else ! double excitation is ab->ab
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(1,1,2)
                     l=exc(1,2,2)
                     rdm2b(i,k,j,l) = rdm2b(i,k,j,l) + cp
                  end if
               end if
            end do
         end do
        end

        subroutine compute_rdm3s(det,Ndet,coef,mo_num,noa,nob,&
                                 Nint,rdm3a,rdm3b,rdm3c,rdm3d)
         integer*8, intent(in)         :: det(Nint,2,Ndet)
         integer, intent(in)           :: Ndet, Nint, mo_num, noa, nob
         double precision, intent(in)  :: coef(Ndet)
         double precision, intent(out) :: rdm3a(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)
         double precision, intent(out) :: rdm3b(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)
         double precision, intent(out) :: rdm3c(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)
         double precision, intent(out) :: rdm3d(mo_num,mo_num,mo_num,mo_num,mo_num,mo_num)

         integer :: i,j,k,l,m,n,o,p,ispin,ishift,ii,jj,kk,ll,mm,nn,oo,pp
         integer :: degree
         integer :: exc(0:3,2,2) ! [(0:n,1,3:mo_idx), hole/part, ispin]
         double precision :: phase, c, cp
         integer :: occ(noa,2)

         rdm3a = 0.d0
         rdm3b = 0.d0
         rdm3c = 0.d0
         rdm3d = 0.d0
         do n=1,Ndet
            call get_occupied(det(:,:,n), Nint, noa, occ)
            c = coef(n)*coef(n)
            !
            ! Diagonal
            !
            ! aaa --->
            do ii=1,noa
               do jj=ii+1,noa
                  do kk=jj+1,noa
                     i = occ(ii,1); j = occ(jj,1); k = occ(kk,1);
                     rdm3a(i,j,k,i,j,k) = rdm3a(i,j,k,i,j,k) + c
                     rdm3a(i,k,j,i,j,k) = rdm3a(i,k,j,i,j,k) - c
                     rdm3a(j,i,k,i,j,k) = rdm3a(j,i,k,i,j,k) - c
                     rdm3a(j,k,i,i,j,k) = rdm3a(j,k,i,i,j,k) + c
                     rdm3a(k,j,i,i,j,k) = rdm3a(k,j,i,i,j,k) - c
                     rdm3a(k,i,j,i,j,k) = rdm3a(k,i,j,i,j,k) + c
                     !
                     rdm3a(i,j,k,i,k,j) = rdm3a(i,j,k,i,k,j) - c
                     rdm3a(i,k,j,i,k,j) = rdm3a(i,k,j,i,k,j) + c
                     rdm3a(j,i,k,i,k,j) = rdm3a(j,i,k,i,k,j) + c
                     rdm3a(j,k,i,i,k,j) = rdm3a(j,k,i,i,k,j) - c
                     rdm3a(k,j,i,i,k,j) = rdm3a(k,j,i,i,k,j) + c
                     rdm3a(k,i,j,i,k,j) = rdm3a(k,i,j,i,k,j) - c
                     !
                     rdm3a(i,j,k,j,i,k) = rdm3a(i,j,k,j,i,k) - c
                     rdm3a(i,k,j,j,i,k) = rdm3a(i,k,j,j,i,k) + c
                     rdm3a(j,i,k,j,i,k) = rdm3a(j,i,k,j,i,k) + c
                     rdm3a(j,k,i,j,i,k) = rdm3a(j,k,i,j,i,k) - c
                     rdm3a(k,j,i,j,i,k) = rdm3a(k,j,i,j,i,k) + c
                     rdm3a(k,i,j,j,i,k) = rdm3a(k,i,j,j,i,k) - c
                     !
                     rdm3a(i,j,k,j,k,i) = rdm3a(i,j,k,j,k,i) + c
                     rdm3a(i,k,j,j,k,i) = rdm3a(i,k,j,j,k,i) - c
                     rdm3a(j,i,k,j,k,i) = rdm3a(j,i,k,j,k,i) - c
                     rdm3a(j,k,i,j,k,i) = rdm3a(j,k,i,j,k,i) + c
                     rdm3a(k,j,i,j,k,i) = rdm3a(k,j,i,j,k,i) - c
                     rdm3a(k,i,j,j,k,i) = rdm3a(k,i,j,j,k,i) + c
                     !
                     rdm3a(i,j,k,k,i,j) = rdm3a(i,j,k,k,i,j) + c
                     rdm3a(i,k,j,k,i,j) = rdm3a(i,k,j,k,i,j) - c
                     rdm3a(j,i,k,k,i,j) = rdm3a(j,i,k,k,i,j) - c
                     rdm3a(j,k,i,k,i,j) = rdm3a(j,k,i,k,i,j) + c
                     rdm3a(k,j,i,k,i,j) = rdm3a(k,j,i,k,i,j) - c
                     rdm3a(k,i,j,k,i,j) = rdm3a(k,i,j,k,i,j) + c
                     !
                     rdm3a(i,j,k,k,j,i) = rdm3a(i,j,k,k,j,i) - c
                     rdm3a(i,k,j,k,j,i) = rdm3a(i,k,j,k,j,i) + c
                     rdm3a(j,i,k,k,j,i) = rdm3a(j,i,k,k,j,i) + c
                     rdm3a(j,k,i,k,j,i) = rdm3a(j,k,i,k,j,i) - c
                     rdm3a(k,j,i,k,j,i) = rdm3a(k,j,i,k,j,i) + c
                     rdm3a(k,i,j,k,j,i) = rdm3a(k,i,j,k,j,i) - c
                  end do
               end do
            end do
            ! aab ---> 
            do ii=1,noa
               do jj=ii+1,noa
                  do kk=1,nob
                     i = occ(ii,1); j = occ(jj,1); k = occ(kk,2);
                     rdm3b(i,j,k,i,j,k) = rdm3b(i,j,k,i,j,k) + c
                     rdm3b(j,i,k,i,j,k) = rdm3b(j,i,k,i,j,k) - c
                     rdm3b(i,j,k,j,i,k) = rdm3b(i,j,k,j,i,k) - c
                     rdm3b(j,i,k,j,i,k) = rdm3b(j,i,k,j,i,k) + c
                  end do
               end do
            end do
            ! abb --->
            do ii=1,noa
               do jj=1,nob
                  do kk=jj+1,nob
                     i = occ(ii,1); j = occ(jj,2); k = occ(kk,2);
                     rdm3c(i,j,k,i,j,k) = rdm3c(i,j,k,i,j,k) + c
                     rdm3c(i,k,j,i,j,k) = rdm3c(i,k,j,i,j,k) - c
                     rdm3c(i,j,k,i,k,j) = rdm3c(i,j,k,i,k,j) - c
                     rdm3c(i,k,j,i,k,j) = rdm3c(i,k,j,i,k,j) + c
                  end do
               end do
            end do
            ! bbb --->
            do ii=1,nob
               do jj=ii+1,nob
                  do kk=jj+1,nob
                     i = occ(ii,2); j = occ(jj,2); k = occ(kk,2);
                     rdm3d(i,j,k,i,j,k) = rdm3d(i,j,k,i,j,k) + c
                     rdm3d(i,k,j,i,j,k) = rdm3d(i,k,j,i,j,k) - c
                     rdm3d(j,i,k,i,j,k) = rdm3d(j,i,k,i,j,k) - c
                     rdm3d(j,k,i,i,j,k) = rdm3d(j,k,i,i,j,k) + c
                     rdm3d(k,j,i,i,j,k) = rdm3d(k,j,i,i,j,k) - c
                     rdm3d(k,i,j,i,j,k) = rdm3d(k,i,j,i,j,k) + c
                     !
                     rdm3d(i,j,k,i,k,j) = rdm3d(i,j,k,i,k,j) - c
                     rdm3d(i,k,j,i,k,j) = rdm3d(i,k,j,i,k,j) + c
                     rdm3d(j,i,k,i,k,j) = rdm3d(j,i,k,i,k,j) + c
                     rdm3d(j,k,i,i,k,j) = rdm3d(j,k,i,i,k,j) - c
                     rdm3d(k,j,i,i,k,j) = rdm3d(k,j,i,i,k,j) + c
                     rdm3d(k,i,j,i,k,j) = rdm3d(k,i,j,i,k,j) - c
                     !
                     rdm3d(i,j,k,j,i,k) = rdm3d(i,j,k,j,i,k) - c
                     rdm3d(i,k,j,j,i,k) = rdm3d(i,k,j,j,i,k) + c
                     rdm3d(j,i,k,j,i,k) = rdm3d(j,i,k,j,i,k) + c
                     rdm3d(j,k,i,j,i,k) = rdm3d(j,k,i,j,i,k) - c
                     rdm3d(k,j,i,j,i,k) = rdm3d(k,j,i,j,i,k) + c
                     rdm3d(k,i,j,j,i,k) = rdm3d(k,i,j,j,i,k) - c
                     !
                     rdm3d(i,j,k,j,k,i) = rdm3d(i,j,k,j,k,i) + c
                     rdm3d(i,k,j,j,k,i) = rdm3d(i,k,j,j,k,i) - c
                     rdm3d(j,i,k,j,k,i) = rdm3d(j,i,k,j,k,i) - c
                     rdm3d(j,k,i,j,k,i) = rdm3d(j,k,i,j,k,i) + c
                     rdm3d(k,j,i,j,k,i) = rdm3d(k,j,i,j,k,i) - c
                     rdm3d(k,i,j,j,k,i) = rdm3d(k,i,j,j,k,i) + c
                     !
                     rdm3d(i,j,k,k,i,j) = rdm3d(i,j,k,k,i,j) + c
                     rdm3d(i,k,j,k,i,j) = rdm3d(i,k,j,k,i,j) - c
                     rdm3d(j,i,k,k,i,j) = rdm3d(j,i,k,k,i,j) - c
                     rdm3d(j,k,i,k,i,j) = rdm3d(j,k,i,k,i,j) + c
                     rdm3d(k,j,i,k,i,j) = rdm3d(k,j,i,k,i,j) - c
                     rdm3d(k,i,j,k,i,j) = rdm3d(k,i,j,k,i,j) + c
                     !
                     rdm3d(i,j,k,k,j,i) = rdm3d(i,j,k,k,j,i) - c
                     rdm3d(i,k,j,k,j,i) = rdm3d(i,k,j,k,j,i) + c
                     rdm3d(j,i,k,k,j,i) = rdm3d(j,i,k,k,j,i) + c
                     rdm3d(j,k,i,k,j,i) = rdm3d(j,k,i,k,j,i) - c
                     rdm3d(k,j,i,k,j,i) = rdm3d(k,j,i,k,j,i) + c
                     rdm3d(k,i,j,k,j,i) = rdm3d(k,i,j,k,j,i) - c
                  end do
               end do
            end do
            !
            ! Off-diagonal
            !
            do m=1,Ndet
               call get_excitation3(det(:,:,n),det(:,:,m),exc,degree,phase,Nint)
               if (degree >= 4 .or. degree==0) cycle
               c = coef(n)*coef(m)
               cp = phase*c
               if (degree==1) then ! single excitation
                  if (exc(0,1,1)==1) then ! single excitation is a->a
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     ! aaa --->
                     do kk=1,noa
                        do ll=kk+1,noa
                           k = occ(kk,1); l = occ(ll,1);
                           rdm3a(j,k,l,i,k,l) = rdm3a(j,k,l,i,k,l) + cp
                           rdm3a(j,k,l,i,l,k) = rdm3a(j,k,l,i,l,k) - cp
                           rdm3a(j,k,l,k,l,i) = rdm3a(j,k,l,k,l,i) + cp
                           rdm3a(j,k,l,k,i,l) = rdm3a(j,k,l,k,i,l) - cp
                           rdm3a(j,k,l,l,i,k) = rdm3a(j,k,l,l,i,k) + cp
                           rdm3a(j,k,l,l,k,i) = rdm3a(j,k,l,l,k,i) - cp
                           !
                           rdm3a(j,l,k,i,k,l) = rdm3a(j,l,k,i,k,l) - cp
                           rdm3a(j,l,k,i,l,k) = rdm3a(j,l,k,i,l,k) + cp
                           rdm3a(j,l,k,k,l,i) = rdm3a(j,l,k,k,l,i) - cp
                           rdm3a(j,l,k,k,i,l) = rdm3a(j,l,k,k,i,l) + cp
                           rdm3a(j,l,k,l,i,k) = rdm3a(j,l,k,l,i,k) - cp
                           rdm3a(j,l,k,l,k,i) = rdm3a(j,l,k,l,k,i) + cp
                           !
                           rdm3a(k,l,j,i,k,l) = rdm3a(k,l,j,i,k,l) + cp
                           rdm3a(k,l,j,i,l,k) = rdm3a(k,l,j,i,l,k) - cp
                           rdm3a(k,l,j,k,l,i) = rdm3a(k,l,j,k,l,i) + cp
                           rdm3a(k,l,j,k,i,l) = rdm3a(k,l,j,k,i,l) - cp
                           rdm3a(k,l,j,l,i,k) = rdm3a(k,l,j,l,i,k) + cp
                           rdm3a(k,l,j,l,k,i) = rdm3a(k,l,j,l,k,i) - cp
                           !
                           rdm3a(k,j,l,i,k,l) = rdm3a(k,j,l,i,k,l) - cp
                           rdm3a(k,j,l,i,l,k) = rdm3a(k,j,l,i,l,k) + cp
                           rdm3a(k,j,l,k,l,i) = rdm3a(k,j,l,k,l,i) - cp
                           rdm3a(k,j,l,k,i,l) = rdm3a(k,j,l,k,i,l) + cp
                           rdm3a(k,j,l,l,i,k) = rdm3a(k,j,l,l,i,k) - cp
                           rdm3a(k,j,l,l,k,i) = rdm3a(k,j,l,l,k,i) + cp
                           !
                           rdm3a(l,j,k,i,k,l) = rdm3a(l,j,k,i,k,l) + cp
                           rdm3a(l,j,k,i,l,k) = rdm3a(l,j,k,i,l,k) - cp
                           rdm3a(l,j,k,k,l,i) = rdm3a(l,j,k,k,l,i) + cp
                           rdm3a(l,j,k,k,i,l) = rdm3a(l,j,k,k,i,l) - cp
                           rdm3a(l,j,k,l,i,k) = rdm3a(l,j,k,l,i,k) + cp
                           rdm3a(l,j,k,l,k,i) = rdm3a(l,j,k,l,k,i) - cp
                           !
                           rdm3a(l,k,j,i,k,l) = rdm3a(l,k,j,i,k,l) - cp
                           rdm3a(l,k,j,i,l,k) = rdm3a(l,k,j,i,l,k) + cp
                           rdm3a(l,k,j,k,l,i) = rdm3a(l,k,j,k,l,i) - cp
                           rdm3a(l,k,j,k,i,l) = rdm3a(l,k,j,k,i,l) + cp
                           rdm3a(l,k,j,l,i,k) = rdm3a(l,k,j,l,i,k) - cp
                           rdm3a(l,k,j,l,k,i) = rdm3a(l,k,j,l,k,i) + cp
                        end do
                     end do
                     ! aab ---> 
                     do kk=1,noa
                        do ll=1,nob
                           k = occ(kk,1); l = occ(ll,2);
                           rdm3b(j,k,l,i,k,l) = rdm3b(j,k,l,i,k,l) + cp
                           rdm3b(j,k,l,k,i,l) = rdm3b(j,k,l,k,i,l) - cp
                           rdm3b(k,j,l,i,k,l) = rdm3b(k,j,l,i,k,l) - cp
                           rdm3b(k,j,l,k,i,l) = rdm3b(k,j,l,k,i,l) + cp
                        end do
                     end do
                     ! abb -->
                     do kk=1,nob
                        do ll=kk+1,nob
                           k = occ(kk,2); l = occ(ll,2);
                           rdm3c(j,k,l,i,k,l) = rdm3c(j,k,l,i,k,l) + cp
                           rdm3c(j,l,k,i,k,l) = rdm3c(j,l,k,i,k,l) - cp
                           rdm3c(j,k,l,i,l,k) = rdm3c(j,k,l,i,l,k) - cp
                           rdm3c(j,l,k,i,l,k) = rdm3c(j,l,k,i,l,k) + cp
                        end do
                     end do
                  else ! single excitation is b->b
                     i=exc(1,1,2)
                     j=exc(1,2,2)
                     ! aab ---> 
                     do kk=1,noa
                        do ll=kk+1,noa
                           k = occ(kk,1); l = occ(ll,1);
                           rdm3b(k,l,j,k,l,i) = rdm3b(k,l,j,k,l,i) + cp
                           rdm3b(k,l,j,l,k,i) = rdm3b(k,l,j,l,k,i) - cp
                           rdm3b(l,k,j,k,l,i) = rdm3b(l,k,j,k,l,i) - cp
                           rdm3b(l,k,j,l,k,i) = rdm3b(l,k,j,l,k,i) + cp
                        end do
                     end do
                     ! abb -->
                     do kk=1,noa
                        do ll=1,nob
                           k = occ(kk,1); l = occ(ll,2);
                           rdm3c(k,j,l,k,i,l) = rdm3c(k,j,l,k,i,l) + cp
                           rdm3c(k,l,j,k,i,l) = rdm3c(k,l,j,k,i,l) - cp
                           rdm3c(k,j,l,k,l,i) = rdm3c(k,j,l,k,l,i) - cp
                           rdm3c(k,l,j,k,l,i) = rdm3c(k,l,j,k,l,i) + cp
                        end do
                     end do
                     ! bbb --->
                     do kk=1,nob
                        do ll=kk+1,nob
                           k = occ(kk,2); l = occ(ll,2);
                           rdm3d(j,k,l,i,k,l) = rdm3d(j,k,l,i,k,l) + cp
                           rdm3d(j,k,l,i,l,k) = rdm3d(j,k,l,i,l,k) - cp
                           rdm3d(j,k,l,k,l,i) = rdm3d(j,k,l,k,l,i) + cp
                           rdm3d(j,k,l,k,i,l) = rdm3d(j,k,l,k,i,l) - cp
                           rdm3d(j,k,l,l,i,k) = rdm3d(j,k,l,l,i,k) + cp
                           rdm3d(j,k,l,l,k,i) = rdm3d(j,k,l,l,k,i) - cp
                           !
                           rdm3d(j,l,k,i,k,l) = rdm3d(j,l,k,i,k,l) - cp
                           rdm3d(j,l,k,i,l,k) = rdm3d(j,l,k,i,l,k) + cp
                           rdm3d(j,l,k,k,l,i) = rdm3d(j,l,k,k,l,i) - cp
                           rdm3d(j,l,k,k,i,l) = rdm3d(j,l,k,k,i,l) + cp
                           rdm3d(j,l,k,l,i,k) = rdm3d(j,l,k,l,i,k) - cp
                           rdm3d(j,l,k,l,k,i) = rdm3d(j,l,k,l,k,i) + cp
                           !
                           rdm3d(k,l,j,i,k,l) = rdm3d(k,l,j,i,k,l) + cp
                           rdm3d(k,l,j,i,l,k) = rdm3d(k,l,j,i,l,k) - cp
                           rdm3d(k,l,j,k,l,i) = rdm3d(k,l,j,k,l,i) + cp
                           rdm3d(k,l,j,k,i,l) = rdm3d(k,l,j,k,i,l) - cp
                           rdm3d(k,l,j,l,i,k) = rdm3d(k,l,j,l,i,k) + cp
                           rdm3d(k,l,j,l,k,i) = rdm3d(k,l,j,l,k,i) - cp
                           !
                           rdm3d(k,j,l,i,k,l) = rdm3d(k,j,l,i,k,l) - cp
                           rdm3d(k,j,l,i,l,k) = rdm3d(k,j,l,i,l,k) + cp
                           rdm3d(k,j,l,k,l,i) = rdm3d(k,j,l,k,l,i) - cp
                           rdm3d(k,j,l,k,i,l) = rdm3d(k,j,l,k,i,l) + cp
                           rdm3d(k,j,l,l,i,k) = rdm3d(k,j,l,l,i,k) - cp
                           rdm3d(k,j,l,l,k,i) = rdm3d(k,j,l,l,k,i) + cp
                           !
                           rdm3d(l,j,k,i,k,l) = rdm3d(l,j,k,i,k,l) + cp
                           rdm3d(l,j,k,i,l,k) = rdm3d(l,j,k,i,l,k) - cp
                           rdm3d(l,j,k,k,l,i) = rdm3d(l,j,k,k,l,i) + cp
                           rdm3d(l,j,k,k,i,l) = rdm3d(l,j,k,k,i,l) - cp
                           rdm3d(l,j,k,l,i,k) = rdm3d(l,j,k,l,i,k) + cp
                           rdm3d(l,j,k,l,k,i) = rdm3d(l,j,k,l,k,i) - cp
                           !
                           rdm3d(l,k,j,i,k,l) = rdm3d(l,k,j,i,k,l) - cp
                           rdm3d(l,k,j,i,l,k) = rdm3d(l,k,j,i,l,k) + cp
                           rdm3d(l,k,j,k,l,i) = rdm3d(l,k,j,k,l,i) - cp
                           rdm3d(l,k,j,k,i,l) = rdm3d(l,k,j,k,i,l) + cp
                           rdm3d(l,k,j,l,i,k) = rdm3d(l,k,j,l,i,k) - cp
                           rdm3d(l,k,j,l,k,i) = rdm3d(l,k,j,l,k,i) + cp
                        end do
                     end do
                  end if
               else if (degree==2) then
                  if (exc(0,1,1)==2) then ! double excitation is aa->aa
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(2,1,1)
                     l=exc(2,2,1)
                     ! aaa --->
                     do oo=1,noa
                        o = occ(oo,1)
                        rdm3a(j,l,o,i,k,o) = rdm3a(j,l,o,i,k,o) + cp
                        rdm3a(j,l,o,i,o,k) = rdm3a(j,l,o,i,o,k) - cp
                        rdm3a(j,l,o,k,o,i) = rdm3a(j,l,o,k,o,i) + cp
                        rdm3a(j,l,o,k,i,o) = rdm3a(j,l,o,k,i,o) - cp
                        rdm3a(j,l,o,o,i,k) = rdm3a(j,l,o,o,i,k) + cp
                        rdm3a(j,l,o,o,k,i) = rdm3a(j,l,o,o,k,i) - cp
                        !
                        rdm3a(j,o,l,i,k,o) = rdm3a(j,o,l,i,k,o) - cp
                        rdm3a(j,o,l,i,o,k) = rdm3a(j,o,l,i,o,k) + cp
                        rdm3a(j,o,l,k,o,i) = rdm3a(j,o,l,k,o,i) - cp
                        rdm3a(j,o,l,k,i,o) = rdm3a(j,o,l,k,i,o) + cp
                        rdm3a(j,o,l,o,i,k) = rdm3a(j,o,l,o,i,k) - cp
                        rdm3a(j,o,l,o,k,i) = rdm3a(j,o,l,o,k,i) + cp
                        !
                        rdm3a(l,o,j,i,k,o) = rdm3a(l,o,j,i,k,o) + cp
                        rdm3a(l,o,j,i,o,k) = rdm3a(l,o,j,i,o,k) - cp
                        rdm3a(l,o,j,k,o,i) = rdm3a(l,o,j,k,o,i) + cp
                        rdm3a(l,o,j,k,i,o) = rdm3a(l,o,j,k,i,o) - cp
                        rdm3a(l,o,j,o,i,k) = rdm3a(l,o,j,o,i,k) + cp
                        rdm3a(l,o,j,o,k,i) = rdm3a(l,o,j,o,k,i) - cp
                        !
                        rdm3a(l,j,o,i,k,o) = rdm3a(l,j,o,i,k,o) - cp
                        rdm3a(l,j,o,i,o,k) = rdm3a(l,j,o,i,o,k) + cp
                        rdm3a(l,j,o,k,o,i) = rdm3a(l,j,o,k,o,i) - cp
                        rdm3a(l,j,o,k,i,o) = rdm3a(l,j,o,k,i,o) + cp
                        rdm3a(l,j,o,o,i,k) = rdm3a(l,j,o,o,i,k) - cp
                        rdm3a(l,j,o,o,k,i) = rdm3a(l,j,o,o,k,i) + cp
                        !
                        rdm3a(o,j,l,i,k,o) = rdm3a(o,j,l,i,k,o) + cp
                        rdm3a(o,j,l,i,o,k) = rdm3a(o,j,l,i,o,k) - cp
                        rdm3a(o,j,l,k,o,i) = rdm3a(o,j,l,k,o,i) + cp
                        rdm3a(o,j,l,k,i,o) = rdm3a(o,j,l,k,i,o) - cp
                        rdm3a(o,j,l,o,i,k) = rdm3a(o,j,l,o,i,k) + cp
                        rdm3a(o,j,l,o,k,i) = rdm3a(o,j,l,o,k,i) - cp
                        !
                        rdm3a(o,l,j,i,k,o) = rdm3a(o,l,j,i,k,o) - cp
                        rdm3a(o,l,j,i,o,k) = rdm3a(o,l,j,i,o,k) + cp
                        rdm3a(o,l,j,k,o,i) = rdm3a(o,l,j,k,o,i) - cp
                        rdm3a(o,l,j,k,i,o) = rdm3a(o,l,j,k,i,o) + cp
                        rdm3a(o,l,j,o,i,k) = rdm3a(o,l,j,o,i,k) - cp
                        rdm3a(o,l,j,o,k,i) = rdm3a(o,l,j,o,k,i) + cp
                     end do
                     ! aab --->
                     do oo=1,nob
                        o = occ(oo,2)
                        rdm3b(j,l,o,i,k,o) = rdm3b(j,l,o,i,k,o) + cp
                        rdm3b(l,j,o,i,k,o) = rdm3b(l,j,o,i,k,o) - cp
                        rdm3b(j,l,o,k,i,o) = rdm3b(j,l,o,k,i,o) - cp
                        rdm3b(l,j,o,k,i,o) = rdm3b(l,j,o,k,i,o) + cp
                     end do
                  else if (exc(0,1,2)==2) then ! double excitation is bb->bb
                     i=exc(1,1,2)
                     j=exc(1,2,2)
                     k=exc(2,1,2)
                     l=exc(2,2,2)
                     ! bbb --->
                     do oo=1,nob
                        o = occ(oo,2)
                        rdm3d(j,l,o,i,k,o) = rdm3d(j,l,o,i,k,o) + cp
                        rdm3d(j,l,o,i,o,k) = rdm3d(j,l,o,i,o,k) - cp
                        rdm3d(j,l,o,k,o,i) = rdm3d(j,l,o,k,o,i) + cp
                        rdm3d(j,l,o,k,i,o) = rdm3d(j,l,o,k,i,o) - cp
                        rdm3d(j,l,o,o,i,k) = rdm3d(j,l,o,o,i,k) + cp
                        rdm3d(j,l,o,o,k,i) = rdm3d(j,l,o,o,k,i) - cp
                        !
                        rdm3d(j,o,l,i,k,o) = rdm3d(j,o,l,i,k,o) - cp
                        rdm3d(j,o,l,i,o,k) = rdm3d(j,o,l,i,o,k) + cp
                        rdm3d(j,o,l,k,o,i) = rdm3d(j,o,l,k,o,i) - cp
                        rdm3d(j,o,l,k,i,o) = rdm3d(j,o,l,k,i,o) + cp
                        rdm3d(j,o,l,o,i,k) = rdm3d(j,o,l,o,i,k) - cp
                        rdm3d(j,o,l,o,k,i) = rdm3d(j,o,l,o,k,i) + cp
                        !
                        rdm3d(l,o,j,i,k,o) = rdm3d(l,o,j,i,k,o) + cp
                        rdm3d(l,o,j,i,o,k) = rdm3d(l,o,j,i,o,k) - cp
                        rdm3d(l,o,j,k,o,i) = rdm3d(l,o,j,k,o,i) + cp
                        rdm3d(l,o,j,k,i,o) = rdm3d(l,o,j,k,i,o) - cp
                        rdm3d(l,o,j,o,i,k) = rdm3d(l,o,j,o,i,k) + cp
                        rdm3d(l,o,j,o,k,i) = rdm3d(l,o,j,o,k,i) - cp
                        !
                        rdm3d(l,j,o,i,k,o) = rdm3d(l,j,o,i,k,o) - cp
                        rdm3d(l,j,o,i,o,k) = rdm3d(l,j,o,i,o,k) + cp
                        rdm3d(l,j,o,k,o,i) = rdm3d(l,j,o,k,o,i) - cp
                        rdm3d(l,j,o,k,i,o) = rdm3d(l,j,o,k,i,o) + cp
                        rdm3d(l,j,o,o,i,k) = rdm3d(l,j,o,o,i,k) - cp
                        rdm3d(l,j,o,o,k,i) = rdm3d(l,j,o,o,k,i) + cp
                        !
                        rdm3d(o,j,l,i,k,o) = rdm3d(o,j,l,i,k,o) + cp
                        rdm3d(o,j,l,i,o,k) = rdm3d(o,j,l,i,o,k) - cp
                        rdm3d(o,j,l,k,o,i) = rdm3d(o,j,l,k,o,i) + cp
                        rdm3d(o,j,l,k,i,o) = rdm3d(o,j,l,k,i,o) - cp
                        rdm3d(o,j,l,o,i,k) = rdm3d(o,j,l,o,i,k) + cp
                        rdm3d(o,j,l,o,k,i) = rdm3d(o,j,l,o,k,i) - cp
                        !
                        rdm3d(o,l,j,i,k,o) = rdm3d(o,l,j,i,k,o) - cp
                        rdm3d(o,l,j,i,o,k) = rdm3d(o,l,j,i,o,k) + cp
                        rdm3d(o,l,j,k,o,i) = rdm3d(o,l,j,k,o,i) - cp
                        rdm3d(o,l,j,k,i,o) = rdm3d(o,l,j,k,i,o) + cp
                        rdm3d(o,l,j,o,i,k) = rdm3d(o,l,j,o,i,k) - cp
                        rdm3d(o,l,j,o,k,i) = rdm3d(o,l,j,o,k,i) + cp
                     end do
                     ! abb --->
                     do oo=1,noa
                        o = occ(oo,1)
                        rdm3c(o,i,k,o,j,l) = rdm3c(o,i,k,o,j,l) + cp
                        rdm3c(o,k,i,o,j,l) = rdm3c(o,k,i,o,j,l) - cp
                        rdm3c(o,i,k,o,l,j) = rdm3c(o,i,k,o,l,j) - cp
                        rdm3c(o,k,i,o,l,j) = rdm3c(o,k,i,o,l,j) + cp
                     end do
                  else ! double excitation is ab->ab
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(1,1,2)
                     l=exc(1,2,2)
                     ! aab --->
                     do oo=1,noa
                        o = occ(oo,1)
                        rdm3b(j,o,l,i,o,k) = rdm3b(j,o,l,i,o,k) + cp
                        rdm3b(j,o,l,o,i,k) = rdm3b(j,o,l,o,i,k) - cp
                        rdm3b(o,j,l,i,o,k) = rdm3b(o,j,l,i,o,k) - cp
                        rdm3b(o,j,l,o,i,k) = rdm3b(o,j,l,o,i,k) + cp
                     end do
                     ! abb --->
                     do oo=1,nob
                        o = occ(oo,2)
                        rdm3c(i,k,o,j,l,o) = rdm3c(i,k,o,j,l,o) + cp
                        rdm3c(i,o,k,j,l,o) = rdm3c(i,o,k,j,l,o) - cp
                        rdm3c(i,k,o,j,o,l) = rdm3c(i,k,o,j,o,l) - cp
                        rdm3c(i,o,k,j,o,l) = rdm3c(i,o,k,j,o,l) + cp
                     end do
                  end if
               else if (degree==3) then
                  if (exc(0,1,1)==3) then ! triple excitation is aaa->aaa
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(2,1,1)
                     l=exc(2,2,1)
                     o=exc(3,1,1)
                     p=exc(3,2,1)
                     ! aaa --->
                     rdm3a(j,l,p,i,k,o) = rdm3a(j,l,p,i,k,o) + cp
                     rdm3a(j,l,p,i,o,k) = rdm3a(j,l,p,i,o,k) - cp
                     rdm3a(j,l,p,k,o,i) = rdm3a(j,l,p,k,o,i) + cp
                     rdm3a(j,l,p,k,i,o) = rdm3a(j,l,p,k,i,o) - cp
                     rdm3a(j,l,p,o,i,k) = rdm3a(j,l,p,o,i,k) + cp
                     rdm3a(j,l,p,o,k,i) = rdm3a(j,l,p,o,k,i) - cp
                     !
                     rdm3a(j,p,l,i,k,o) = rdm3a(j,p,l,i,k,o) - cp
                     rdm3a(j,p,l,i,o,k) = rdm3a(j,p,l,i,o,k) + cp
                     rdm3a(j,p,l,k,o,i) = rdm3a(j,p,l,k,o,i) - cp
                     rdm3a(j,p,l,k,i,o) = rdm3a(j,p,l,k,i,o) + cp
                     rdm3a(j,p,l,o,i,k) = rdm3a(j,p,l,o,i,k) - cp
                     rdm3a(j,p,l,o,k,i) = rdm3a(j,p,l,o,k,i) + cp
                     !
                     rdm3a(l,p,j,i,k,o) = rdm3a(l,p,j,i,k,o) + cp
                     rdm3a(l,p,j,i,o,k) = rdm3a(l,p,j,i,o,k) - cp
                     rdm3a(l,p,j,k,o,i) = rdm3a(l,p,j,k,o,i) + cp
                     rdm3a(l,p,j,k,i,o) = rdm3a(l,p,j,k,i,o) - cp
                     rdm3a(l,p,j,o,i,k) = rdm3a(l,p,j,o,i,k) + cp
                     rdm3a(l,p,j,o,k,i) = rdm3a(l,p,j,o,k,i) - cp
                     !
                     rdm3a(l,j,p,i,k,o) = rdm3a(l,j,p,i,k,o) - cp
                     rdm3a(l,j,p,i,o,k) = rdm3a(l,j,p,i,o,k) + cp
                     rdm3a(l,j,p,k,o,i) = rdm3a(l,j,p,k,o,i) - cp
                     rdm3a(l,j,p,k,i,o) = rdm3a(l,j,p,k,i,o) + cp
                     rdm3a(l,j,p,o,i,k) = rdm3a(l,j,p,o,i,k) - cp
                     rdm3a(l,j,p,o,k,i) = rdm3a(l,j,p,o,k,i) + cp
                     !
                     rdm3a(p,j,l,i,k,o) = rdm3a(p,j,l,i,k,o) + cp
                     rdm3a(p,j,l,i,o,k) = rdm3a(p,j,l,i,o,k) - cp
                     rdm3a(p,j,l,k,o,i) = rdm3a(p,j,l,k,o,i) + cp
                     rdm3a(p,j,l,k,i,o) = rdm3a(p,j,l,k,i,o) - cp
                     rdm3a(p,j,l,o,i,k) = rdm3a(p,j,l,o,i,k) + cp
                     rdm3a(p,j,l,o,k,i) = rdm3a(p,j,l,o,k,i) - cp
                     !
                     rdm3a(p,l,j,i,k,o) = rdm3a(p,l,j,i,k,o) - cp
                     rdm3a(p,l,j,i,o,k) = rdm3a(p,l,j,i,o,k) + cp
                     rdm3a(p,l,j,k,o,i) = rdm3a(p,l,j,k,o,i) - cp
                     rdm3a(p,l,j,k,i,o) = rdm3a(p,l,j,k,i,o) + cp
                     rdm3a(p,l,j,o,i,k) = rdm3a(p,l,j,o,i,k) - cp
                     rdm3a(p,l,j,o,k,i) = rdm3a(p,l,j,o,k,i) + cp
                  else if (exc(0,1,1)==2) then ! triple excitation is abb->abb
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(2,1,1)
                     l=exc(2,2,1)
                     o=exc(1,1,2)
                     p=exc(1,2,2)
                     ! aab --->
                     rdm3b(j,l,p,i,k,o) = rdm3b(j,l,p,i,k,o) + cp
                     rdm3b(j,l,p,k,i,o) = rdm3b(j,l,p,k,i,o) - cp
                     rdm3b(l,j,p,i,k,o) = rdm3b(l,j,p,i,k,o) - cp
                     rdm3b(l,j,p,k,i,o) = rdm3b(l,j,p,k,i,o) + cp
                  else if (exc(0,1,1)==1) then ! triple excitation is abb->abb
                     i=exc(1,1,1)
                     j=exc(1,2,1)
                     k=exc(1,1,2)
                     l=exc(1,2,2)
                     o=exc(2,1,2)
                     p=exc(2,2,2)
                     ! abb --->
                     rdm3c(j,l,p,i,k,o) = rdm3c(j,l,p,i,k,o) + cp
                     rdm3c(j,l,p,i,o,k) = rdm3c(j,l,p,i,o,k) - cp
                     rdm3c(j,p,l,i,k,o) = rdm3c(j,p,l,i,k,o) - cp
                     rdm3c(j,p,l,i,o,k) = rdm3c(j,p,l,i,o,k) + cp
                  else ! triple excitation is bbb->bbb
                     i=exc(1,1,2)
                     j=exc(1,2,2)
                     k=exc(2,1,2)
                     l=exc(2,2,2)
                     o=exc(3,1,2)
                     p=exc(3,2,2)
                     ! bbb --->
                     rdm3d(j,l,p,i,k,o) = rdm3d(j,l,p,i,k,o) + cp
                     rdm3d(j,l,p,i,o,k) = rdm3d(j,l,p,i,o,k) - cp
                     rdm3d(j,l,p,k,o,i) = rdm3d(j,l,p,k,o,i) + cp
                     rdm3d(j,l,p,k,i,o) = rdm3d(j,l,p,k,i,o) - cp
                     rdm3d(j,l,p,o,i,k) = rdm3d(j,l,p,o,i,k) + cp
                     rdm3d(j,l,p,o,k,i) = rdm3d(j,l,p,o,k,i) - cp
                     !
                     rdm3d(j,p,l,i,k,o) = rdm3d(j,p,l,i,k,o) - cp
                     rdm3d(j,p,l,i,o,k) = rdm3d(j,p,l,i,o,k) + cp
                     rdm3d(j,p,l,k,o,i) = rdm3d(j,p,l,k,o,i) - cp
                     rdm3d(j,p,l,k,i,o) = rdm3d(j,p,l,k,i,o) + cp
                     rdm3d(j,p,l,o,i,k) = rdm3d(j,p,l,o,i,k) - cp
                     rdm3d(j,p,l,o,k,i) = rdm3d(j,p,l,o,k,i) + cp
                     !
                     rdm3d(l,p,j,i,k,o) = rdm3d(l,p,j,i,k,o) + cp
                     rdm3d(l,p,j,i,o,k) = rdm3d(l,p,j,i,o,k) - cp
                     rdm3d(l,p,j,k,o,i) = rdm3d(l,p,j,k,o,i) + cp
                     rdm3d(l,p,j,k,i,o) = rdm3d(l,p,j,k,i,o) - cp
                     rdm3d(l,p,j,o,i,k) = rdm3d(l,p,j,o,i,k) + cp
                     rdm3d(l,p,j,o,k,i) = rdm3d(l,p,j,o,k,i) - cp
                     !
                     rdm3d(l,j,p,i,k,o) = rdm3d(l,j,p,i,k,o) - cp
                     rdm3d(l,j,p,i,o,k) = rdm3d(l,j,p,i,o,k) + cp
                     rdm3d(l,j,p,k,o,i) = rdm3d(l,j,p,k,o,i) - cp
                     rdm3d(l,j,p,k,i,o) = rdm3d(l,j,p,k,i,o) + cp
                     rdm3d(l,j,p,o,i,k) = rdm3d(l,j,p,o,i,k) - cp
                     rdm3d(l,j,p,o,k,i) = rdm3d(l,j,p,o,k,i) + cp
                     !
                     rdm3d(p,j,l,i,k,o) = rdm3d(p,j,l,i,k,o) + cp
                     rdm3d(p,j,l,i,o,k) = rdm3d(p,j,l,i,o,k) - cp
                     rdm3d(p,j,l,k,o,i) = rdm3d(p,j,l,k,o,i) + cp
                     rdm3d(p,j,l,k,i,o) = rdm3d(p,j,l,k,i,o) - cp
                     rdm3d(p,j,l,o,i,k) = rdm3d(p,j,l,o,i,k) + cp
                     rdm3d(p,j,l,o,k,i) = rdm3d(p,j,l,o,k,i) - cp
                     !
                     rdm3d(p,l,j,i,k,o) = rdm3d(p,l,j,i,k,o) - cp
                     rdm3d(p,l,j,i,o,k) = rdm3d(p,l,j,i,o,k) + cp
                     rdm3d(p,l,j,k,o,i) = rdm3d(p,l,j,k,o,i) - cp
                     rdm3d(p,l,j,k,i,o) = rdm3d(p,l,j,k,i,o) + cp
                     rdm3d(p,l,j,o,i,k) = rdm3d(p,l,j,o,i,k) - cp
                     rdm3d(p,l,j,o,k,i) = rdm3d(p,l,j,o,k,i) + cp
                  end if
               end if
            end do
         end do
        end

        
        subroutine sort_determinants(sigma,coef,det,N_int,ndet,A,nloc,ispin,e_diag)

              integer, intent(in) :: N_int, ndet, nloc, ispin

              integer, intent(inout) :: A(nloc+1)
              integer*8, intent(inout) :: det(N_int,2,ndet)
              double precision, intent(inout) :: coef(ndet)
              double precision, intent(inout) :: sigma(ndet)
              double precision, intent(inout) :: e_diag(ndet)

              integer :: j, p
              integer :: idet, iloc, p1, p2
              integer, allocatable :: idx(:)
              double precision, allocatable :: temp(:)
              
              !!! ASSUMING N_int=1 !!!
              j = 1

              ! get the sorting array for ispin-major sorting (alpha or beta)
              allocate(idx(ndet),temp(ndet))
!              do j=1,N_int
!                 call argsort_i8(det(j,ispin,:), idx)
!                 ! apply sorting array to determinants, CI coefficients, and sigma array
!                 !det(j,:,:) = det(j,:,idx)
!!                 temp(:) = coef(:)
!!                 do p=1,ndet
!!                    coef(p) = temp(idx(p))
!!                 end do
!!                 temp(:) = sigma(:)
!!                 do p=1,ndet
!!                    sigma(p) = temp(idx(p))
!!                 end do
!!                 temp(:) = e_diag(:)
!!                 do p=1,ndet
!!                    e_diag(p) = temp(idx(p))
!!                 end do
!              end do
              call argsort_i8(det(j,ispin,:), idx)
              ! apply sorting array to determinants, CI coefficients, and sigma array
              det(j,:,:) = det(j,:,idx)
              temp(:) = coef(:)
              do p=1,ndet
                 coef(p) = temp(idx(p))
              end do
              temp(:) = sigma(:)
              do p=1,ndet
                 sigma(p) = temp(idx(p))
              end do
              temp(:) = e_diag(:)
              do p=1,ndet
                 e_diag(p) = temp(idx(p))
              end do
              deallocate(idx,temp)
              ! create array A, where A(n) marks the first occurence of a new spin bucket
              iloc = 1
              A(1) = 1
              A(nloc+1) = ndet + 1
              do idet=1,ndet-1
                 ! get consecutive lexcial indices
                 p1 = det(j,ispin,idet)
                 p2 = det(j,ispin,idet+1)
                 ! if change occurs between consecutive indices, record these locations in loc_arr as new start/end points
                 if (p1 /= p2) then
                    iloc = iloc+1
                    A(iloc) = idet+1
                 end if
              end do
      end subroutine sort_determinants
      
      subroutine sort_determinants_simple(det,H,N_int,ndet,A,nloc,ispin)

              integer, intent(in) :: N_int, ndet, nloc, ispin

              integer, intent(inout) :: A(nloc+1)
              integer*8, intent(inout) :: det(N_int,2,ndet)
              double precision, intent(inout) :: H(ndet,ndet)

              integer :: j, p, q
              integer :: idet, iloc, p1, p2
              integer, allocatable :: idx(:)
              double precision, allocatable :: htemp(:,:)
              
              !!! ASSUMING N_int=1 !!!
              j = 1

              ! get the sorting array for ispin-major sorting (alpha or beta)
              allocate(idx(ndet),htemp(ndet,ndet))
              call argsort_i8(det(j,ispin,:), idx)
              ! apply sorting array to determinants, CI coefficients, and sigma array
              det(j,:,:) = det(j,:,idx)
              ! store old H matrix
              htemp(:,:) = H(:,:)
              do p=1,ndet
                 do q=1,ndet
                    !print*, 'hnew(',p,q,')=','hold(',idx(p),idx(q),')'
                    H(p,q)=htemp(idx(p),idx(q))
                    !H(idx(p),idx(q)) = H(p,q)
                 end do
              end do
              deallocate(idx,htemp)
              ! create array A, where A(n) marks the first occurence of a new spin bucket
              iloc = 1
              A(1) = 1
              A(nloc+1) = ndet + 1
              do idet=1,ndet-1
                 ! get consecutive lexcial indices
                 p1 = det(j,ispin,idet)
                 p2 = det(j,ispin,idet+1)
                 ! if change occurs between consecutive indices, record these locations in loc_arr as new start/end points
                 if (p1 /= p2) then
                    iloc = iloc+1
                    A(iloc) = idet+1
                 end if
              end do
      end subroutine sort_determinants_simple

      subroutine argsort_i8(r,d)

              integer*8, intent(in), dimension(:) :: r
              integer, intent(out), dimension(size(r)) :: d

              integer, dimension(size(r)) :: il

              integer :: stepsize
              integer :: i, j, n, left, k, ksize

              n = size(r)

              do i=1,n
                 d(i)=i
              end do

              if (n==1) return

              stepsize = 1
              do while (stepsize < n)
                 do left = 1, n-stepsize,stepsize*2
                    i = left
                    j = left+stepsize
                    ksize = min(stepsize*2,n-left+1)
                    k=1

                    do while (i < left+stepsize .and. j < left+ksize)
                       if (r(d(i)) < r(d(j))) then
                          il(k) = d(i)
                          i = i+1
                          k = k+1
                       else
                          il(k) = d(j)
                          j = j+1
                          k = k+1
                       endif
                    enddo

                    if (i < left+stepsize) then
                       ! fill up remaining from left
                       il(k:ksize) = d(i:left+stepsize-1)
                    else
                       ! fill up remaining from right
                       il(k:ksize) = d(j:left+ksize-1)
                    endif
                    d(left:left+ksize-1) = il(1:ksize)
                 end do
                 stepsize = stepsize*2
              end do

      end subroutine argsort_i8
   
      subroutine check_spin_buckets( N_int, ndet, nloc, A, det )
         
            integer, intent(in) :: N_int, ndet, nloc
            integer*8, intent(in) :: det(N_int,2,ndet)
            integer, intent(in) :: A(nloc+1)
            
            integer :: b1, b2
            
            do b1=1,nloc
               print*, 'spin bucket', b1
               do b2=A(b1),A(b1+1)-1
                  print*, det(1,:,b2)
               end do
            end do
         
      end subroutine check_spin_buckets
   
end module ci
