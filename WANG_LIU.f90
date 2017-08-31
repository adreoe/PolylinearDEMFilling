! ��������� Wang_Liu � ���������������� ������� ��������� �����-�� ��� ���������� ���������
! ��� ��������� ������������ �������� ��� (Z), ��������� � ������������ ������� (Z_flat), � ����� ������� � ����������� ����������� � �������� ��������� (depr)
    
    subroutine Wang_Liu(Z,Z_flat,depr,out_list,q_out,Zmax,Nx,Ny,NODATA)
    
    implicit none
    
    ! ����������, ������������ � ������������
    real :: Z(Nx,Ny) ! �������� ���
    real :: Z_flat(Nx,Ny) ! ��������� ���������� (��� ���������)
    integer :: depr(Nx,Ny) ! ������������� ���������� �� ������. 1 - ���������, 0 - �� ���������.
    integer :: out_list(Nx*Ny) ! ������ ����� ������ (����������, ���������� �� ����������� ������)
    integer :: q_out ! �������� ��� ������ ����� ������
    real :: Zmax ! ������������ ������ � �������� ������
    integer :: Nx,Ny ! ����������� �����
    real :: NODATA ! �������� ���� �������

    ! ���������� ���������� ������������
    integer, allocatable :: mask(:,:) ! ������, � ������� ����� ���������� ���������������/����������������� ������. 0 - �� �����������, 1 - � �������, 2 - �����������
    integer, allocatable :: list_cells(:) ! ����� ������ �����
    real, allocatable :: list_Z(:) ! ������ ����� �����. ���������� ������ list_cells. ��������� ��� �������� ��������� ������ (�� ���� ������ ��� ���������� � �������� �������)
    integer :: i,j ! ��������� � ������
    integer :: q1 ! ����� ���������� ��������� � ������ �����
    integer :: q2, q3 ! ��������� � ������ (������� ����� � ����� ������)
    integer :: q_max ! ����� ���������� ����� ��� ���������. � ������ ������ ��������� ����� ���������� ����� � �������, ����� ��������� (�� ���� �������� "��� ������")
    integer :: index
    real :: Ztres
    integer :: percent_complete, percent_0
    integer :: c1,c2,r1,r2
    
    !! ���������� ���� 3�3 (������� �� ����������� "������" ������ ������� �������)
    !integer :: k
    !integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    !integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    
    ! ���������� ���� 3�3 (������� �������� ������ �������, ����� ������� �� ���������)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! ������ ������
    ! ���������� ����������� ���������� � ��������.
    q_max = Nx*Ny
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny)
    
    depr(1:Nx,1:Ny) = 0
    allocate(mask(Nx,Ny))
    mask(1:Nx,1:Ny) = 0
    allocate(list_cells(Nx*Ny))
    allocate(list_Z(Nx*Ny))
    
    percent_0 = 0
    
    ! ����������� ��������������� ������ �����
    q2 = 0
    call initialize_set()
    
    ! ��������� ���
    q2 = 1
    q_out = 0
    do 
        if (q2 > q_max) exit
        call one_to_two(c1,r1,list_cells(q2)) ! �������� ������ ��� ��������� �� ������� ����
    
        do k = 1,8 ! ��� � ����� � �������
            c2 = c1+kx(k); r2 = r1+ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (mask(c2,r2) /= 0) cycle ! ���� � ����� ������ �� ��� ����, ��������� �� ��������� ��������
            
            ! ��������� �������������� ������ � ������
            call add_cell(c2,r2,q2)
                       
            ! ���������, ���� �� � ������ ���������
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
            Z_flat(c2,r2) = Z_flat(c1,r1)
            depr(c2,r2) = 1
            
            if( depr(c1,r1) == 0) then ! ���� �����, �� ������� �� ������, �� �������� ����������, ��������� � ��� ����� ������
                depr(c1,r1) = -2
                q_out = q_out + 1
                out_list(q_out) = two_to_one(c1,r1)
            endif
            
            endif
        enddo
        mask(c1,r1) = 2
        q2 = q2+1
        ! ����������� �������� ����������
        percent_complete = int(real(q2)/real(q_max) * 100)
        if (percent_complete > percent_0) then
            print *, "Searching for pits:", percent_complete, "% completed"
            percent_0 = percent_complete
        endif 
    enddo
    
    print *, "Found", q_out, "depressions"
    if(allocated(mask)) deallocate(mask)
    if(allocated(list_cells)) deallocate(list_cells)
    if(allocated(list_Z)) deallocate(list_Z)
    
    contains

! ������������ ������������ ��������������� ������
subroutine initialize_set()
! ���� �����:
! ��������������� ��� ������ ���. �� ��������� ���������, ��� ������ �� ������� ��������� � �������������� ������ (decision = .false.)
! ������, ���� ������ ��������� �� �������, ��� � �� ����� "��� ������" � ������� �������� �� ������������� (decision = .true.)
! ������ �� ���������� "��� ������" � ������ �� ���������, � ����� ��������� �� ������������
! � ����� ������������, ���� ������� �������������, ������ ����������� � �������������� ������

logical :: decision

q1 = 0
! ������������� �� ������� ��� ������ �������
do c1 = 1,Nx
    do r1 = 1,Ny
        decision = .false.
        ! ���� � ������ ��� ������, ��������� � �� ������������ ��������
        if (Z(c1,r1) == nodata) then
            mask(c1,r1) = 2
            q_max = q_max - 1
            cycle
        endif
        ! ���� � ������ ���� ������, ��������������� � ��� ��������������
        do k = 1,8
            c2 = c1+kx(k); r2 = r1+ky(k) ! ���� � ������
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) then ! ���� ����� �� ���������� �����, ��������� ������� ������� ������ � ������
                decision = .true.
                exit
            endif
            if (Z(c2,r2) == nodata) then ! ���� ����� ����������, �� � �� "��� ������", ��������� ������� ������� ������ � ������
                decision = .true.
            endif
        enddo
        ! ���������� ������ � ������
        if (decision == .true.) then
            call add_cell(c1,r1,0)
        endif     
        
    enddo
enddo


end subroutine initialize_set

! ������������ ���������� �������� � ������
subroutine add_cell(c,r,q2)

integer :: c,r,q2 ! ������� � ������ ������, ����������� � ������, � ����� ����� �������� �������� � ������
real(4) :: Z_index
integer :: index
integer :: q
integer :: new_element ! ������, ��� ������� ����� ������� ����� �������. �� ��������� � � ����� ������
integer :: first, last
integer :: ncase
    
    index = two_to_one(c,r)

    q1 = q1+1
    new_element = q1 ! ����������, ��� ���, "�� ���������"
    
    if (q1 > 1) then
        if (Z(c,r) < List_z(q1-1)) then ! ���� ����� ������� ������ ����������, ��������� ��������� ������.
            ! ���� ����� ������� ������ ����������, ����� ������ ��� � �����
            first = q2+1; last = q1-1
            
            ! ������� �����������
            if (last /= first) then
                do while (last /= first)
                    if (Z(c,r) < List_z( (last-first)/2 + first) ) then
                        last = (last-first)/2 + first
                    else
                        first = (last-first)/2 + first + 1
                    endif
                enddo
            endif
            new_element = first
        call move_cells_in_queue(new_element,q1)
        endif
        
    endif        
    
    ! ���������� ����� ������� � ��� ������ � �� �������, ������� ���� ���������� (q, ��� �� new_element)
    mask(c,r) = 1
    list_cells(new_element) = index
    list_Z(new_element) = Z(c,r)
end subroutine add_cell

! ������������ ����������� ��������� � ������
subroutine move_cells_in_queue(first_cell,last_cell)
integer :: number
integer :: first_cell,last_cell

do number = last_cell, first_cell+1, -1
    list_cells(number) = list_cells(number-1)
    list_Z(number) = list_Z(number-1)
enddo

end subroutine move_cells_in_queue

! ������� ��� ����������� ���������� ������� � ����������
integer function two_to_one(column,row)
integer :: row,column
two_to_one = (column-1)*Ny+row
end function two_to_one

! ������������ ��� �������������� ���������� ������� �� �����������
subroutine one_to_two(column,row,index)
integer :: column,row,index
row = mod(index,Ny)
if (row == 0) then ! ���� �� ������ ����� ��������, ������ ����� ������������ �����������
    row = Ny
    index = index - Ny
endif
column = index/Ny + 1    
end subroutine one_to_two

end subroutine Wang_Liu