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
                        do idet=1,ndet
                           call get_occupied(det(:,:,idet), N_int, noa, occ)
                           do jdet=1,ndet
                              hmatel = 0.0d0
                              call get_excitation(det(:,:,idet),det(:,:,jdet),exc,degree,phase,N_int)
                              if (degree==0) then
                                 call slater0(occ,noa,nob,mo_num,e1int,e2int,hmatel)
                              else if (degree==1) then
                                 call slater1(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                              else if (degree==2) then
                                 call slater2(occ,noa,nob,mo_num,e1int,e2int,exc,phase,hmatel)
                              end if
                              sigma(idet) = sigma(idet)+hmatel*coef(jdet)
                           end do
                           !print*, "idet = ", idet, "sigma = ", sigma(idet)
                        end do
                end subroutine calc_sigma


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
                 implicit none
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

end module ci
