subroutine Wang_Liu2(Z,Z_flat,depr,out_list,q_out,Zmax,Nx,Ny,NODATA)
    
    use priority_queue_mod
    
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

    !! ���������� ���������� ������������
    integer, allocatable :: mask(:,:) ! ������, � ������� ����� ���������� ���������������/����������������� ������. 
    ! 0 - �� �����������, 1 - � �������, 2 - �����������
    integer :: q2, q_max ! ����� ���������� ����� ��� ���������. 
    ! � ������ ������ ��������� ����� ���������� ����� � �������, ����� ��������� (�� ���� �������� "��� ������")
    integer :: percent_complete, percent_0
    integer :: c1,c2,r1,r2
    real :: Z1,Z2
    
    !! ���������� ���� 3�3 (������� �� ����������� "������" ������ ������� �������)
    !integer :: k
    !integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    !integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    
    ! ���������� ���� 3�3 (������� �������� ������ �������, ����� ������� �� ���������)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
    type (queue) :: q
    type (node)  :: x
    
! ������ ������
    ! ���������� ����������� ���������� � ��������.
    q_max = Nx*Ny
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny)
    
    depr(1:Nx,1:Ny) = 0
    allocate(mask(Nx,Ny))
    mask(1:Nx,1:Ny) = 0
    
    percent_0 = 0
    
    ! ����������� ��������������� ������ �����
    q2 = 0
    call initialize_set()
    
    ! ��������� ���
    q_out = 0
    do while (q%n >0)
        ! ��������� ������ ������� �� ������
        x = q%top()
        Z1 = x%priority; c1 = x%c; r1 = x%r
        do k = 1,8
            c2 = c1+kx(k); r2 = r1 + ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle ! ���� ����� �� ������� ������, ������������
            if (mask(c2,r2) /= 0) cycle ! ������������ ������, � ������� �������� mask �� ����� ���� (������������� �����, Nodata, etc.)
            call q%enqueue(Z(c2,r2), c2, r2)! ��������� � �������
            mask(c2,r2) = 1 ! ������������� �����
            
            ! ���������, ���� �� � ������ ���������
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
                Z_flat(c2,r2) = Z_flat(c1,r1)
                depr(c2,r2) = 1
                ! ���� �����, �� ������� �� ������, �� �������� ����������, ��������� � ��� ����� ������
                if( depr(c1,r1) == 0) then 
                    depr(c1,r1) = -2
                    q_out = q_out + 1
                    out_list(q_out) = (c1 - 1) * Ny + r1
                endif
            endif
        enddo
            
        mask(c1,r1) = 2
            
        ! ����������� �������� ����������
        q2 = q2 + 1
        percent_complete = int(real(q2)/real(q_max) * 100)
        if (percent_complete > percent_0) then
            print *, "Searching for pits:", percent_complete, "% completed"
            percent_0 = percent_complete
        endif 
            
        
            
    !    if (q2 > q_max) exit
    !    call one_to_two(c1,r1,list_cells(q2)) ! �������� ������ ��� ��������� �� ������� ����
    !
    !    do k = 1,8 ! ��� � ����� � �������
    !        c2 = c1+kx(k); r2 = r1+ky(k)
    !        if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
    !        if (mask(c2,r2) /= 0) cycle ! ���� � ����� ������ �� ��� ����, ��������� �� ��������� ��������
    !        
    !        ! ��������� �������������� ������ � ������
    !        call add_cell(c2,r2,q2)
    !                   
    !        ! ���������, ���� �� � ������ ���������
    !        if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
    !        Z_flat(c2,r2) = Z_flat(c1,r1)
    !        depr(c2,r2) = 1
    !        
    !        if( depr(c1,r1) == 0) then ! ���� �����, �� ������� �� ������, �� �������� ����������, ��������� � ��� ����� ������
    !            depr(c1,r1) = -2
    !            q_out = q_out + 1
    !            out_list(q_out) = two_to_one(c1,r1)
    !        endif
    !        
    !        endif
    !    enddo
    !    mask(c1,r1) = 2
    !    q2 = q2+1
    !    ! ����������� �������� ����������
    !    percent_complete = int(real(q2)/real(q_max) * 100)
    !    if (percent_complete > percent_0) then
    !        print *, "Searching for pits:", percent_complete, "% completed"
    !        percent_0 = percent_complete
    !    endif 
    enddo
    if(allocated(mask)) deallocate(mask)

contains

subroutine initialize_set()
! ���� �����:
! ��������������� ��� ������ ���. �� ��������� ���������, ��� ������ �� ������� ��������� � �������������� ������ (decision = .false.)
! ������, ���� ������ ��������� �� �������, ��� � �� ����� "��� ������" � ������� �������� �� ������������� (decision = .true.)
! ������ �� ���������� "��� ������" � ������ �� ���������, � ����� ��������� �� ������������
! � ����� ������������, ���� ������� �������������, ������ ����������� � �������������� ������

logical :: decision
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
            call q%enqueue(Z(c1,r1), c1, r1)
        endif     
    enddo
enddo
end subroutine initialize_set

end subroutine Wang_Liu2