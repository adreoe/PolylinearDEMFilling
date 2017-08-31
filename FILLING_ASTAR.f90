!!!!-----------------------------------!!!!
! ��������� ������������� ��� ���������� ��������� ��������� ����� ������������ ���������� ��������� �*
    
subroutine filling_ASTAR(input_model,depr,col_out,row_out,output_model,Nx,Ny,Zmax,StepX,StepY,nodata)
    implicit none
    
    ! ����������, ������������ � �����������
    real :: input_model(Nx,Ny) ! �������� ������� �����
    real :: output_model(Nx,Ny) ! ������� ��� ������ ���������� ����������
    integer :: depr(Nx,Ny) ! ������� ���������, ������ � ����� ������
    integer :: Nx,Ny ! ����������� ������
    real(8) :: StepX,StepY ! ��� ����� (� �������������, ��� ���� �� X � Y ����� �����������)
    real :: Zmax ! ������������ ������ � �������� �������� ������. ������ ���� ��������� � ��������� �����
    integer :: row_out, col_out ! ������� � ������� ������ ������
    real :: nodata
    logical :: connect_depressions
    
    ! ����������, ���������� � ���������
    real, allocatable :: waylength(:,:) ! ������, � ������� ������������ ����� ����. ��� �� �� �� �����, ��� ��������� ����
    integer, allocatable :: way_prev(:,:) ! � ���� ������ �� �����, �� ����� ����� �� ������������
    integer, allocatable  :: flag_status(:,:) ! ������, ���������� ���������� ������: 
    ! 0 � ��� ��������, 
    ! 1 � ��������� � ���������, ����������, 
    ! 2 � ��������� � �������, ���������� � ������� �� ��������� ���,
    ! 3 � ��������� � ���������, ��������� � ������� �� ���������
    ! -1 � ����� �� ���������
    integer, allocatable :: list_cells(:) ! ����� ������ �����
    real, allocatable :: list_Z(:) ! ������ ����� �����. ! ���������� ������ list_cells. ��������� ��� �������� ��������� ������ (�� ���� ������ ��� ���������� � �������� �������)
    integer, allocatable :: list_borders(:) ! ����� ������ ��������� �����. �����, ����� �� ����������� ������ ��� ��� ������� ��� ���������� �������� �����
    
    integer :: c1,c2,c3,r1,r2,r3 ! ���������� ��� �������� � �������
    real :: length
    real :: dZ ! �������� ����� (������������ ��� ������� ����� ������ �����
    
    ! ������ ���������� ������������� �������, ������ ��������� "���������������" (1 - � �������, 0 - �����������)
    integer :: q1, q2, q3 ! ������� ������ ���������� ����� � ������ ������� ����� � ������
    integer :: s,s2 ! ���������� ��� ����������� ������-������� � ���� ����� ����� (� ����������� ������ ����� ����� � ������)
    integer :: i,j,k ! ��������� ��� ������
    integer :: i1, j1
    
    ! ���������� ��� ������ � ���������� A*
    type (queue) :: queue_points, queue_borders
    type (node)  :: x
    
    
    real(8) :: Z0,Z1,Z2,Z3,Z_temp ! ������ ����� ���������, ��������� ����� ������ ����������� �����
    
    ! ���������� ��� ���������� ������ ��� �������, ����� ������� �� �������� ���������� ����
    real :: dist1, dist2
    
    ! ���������� ���� 3�3 (������� �������� ������ �������, ����� ������� �� ���������)
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! ������ ������
    
    allocate(waylength(Nx,Ny))
    allocate(way_prev(Nx,Ny))
    allocate(flag_status(Nx,Ny))
    allocate(list_cells(Nx*Ny))
    allocate(list_Z(Nx*Ny))
    allocate(list_borders(Nx*Ny))
    
    connect_depressions = .TRUE.
    
    ! ��������� �������� ������
    do i = 1,Nx
        do j = 1,Ny
            if (depr(i,j) > 0) then
                output_model(i,j) = nodata
            else
                output_model(i,j) = input_model(i,j)
            endif
        enddo
    enddo
    
    ! ���������� ������� ������-������ � ����������, ������� ����� ������������ � �����
    c1 = col_out; r1 = row_out; 
    
    ! ���������������� ��������
    list_cells(1:Nx*Ny) = 0 ! ����� ������ �����
    list_Z(1:Nx*Ny) = nodata ! ����� ������ ����� �����
    list_borders(1:Nx*Ny) = 0 ! ����� ������ ����� ������
    waylength(1:Nx,1:Ny) = nodata

    ! ������ ������ ����� � �������
    flag_status(1:Nx,1:Ny) = 0
    q1 = 1; q2 = 1; q3 = 0
    list_cells(1) = two_to_one(c1,r1); list_Z(1) = input_model(c1,r1); flag_status(c1,r1) = -1
    ! ���������� Z1 ������ ������ ������, ������������ � ������������
    Z1 = Zmax ! ���������� �������� ������ ����� ������� � �������� ������ ���������
    waylength(c1,r1) = 0.0
    ! ���������������� ���� �� ���� ��������
    
    ! ��������� ������ ���������� ����� �� ��������� AT
    do
        if (q2 > 1 .and. list_cells(q2) == 0) exit ! ��������� ����� ���������� �����, ����� � ������� �� ������� �����
        ! �� ���� �� ��������� ����� � ������-������� � ������ ������
                
        call one_to_two(c1,r1,list_cells(q2)) ! �������� ������ ��� ��������� �� ������� ����
        !!!x = queue_points%top()
        
        if (depr(c1,r1) == 0 .or. depr(c1,r1) == -1) cycle ! ���������� ������ �� �������, ���� ��� �� ��������� � ���������
        
        do k = 1,8 ! ��� � ����� � �������
            c2 = c1+kx(k); r2 = r1+ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (depr(c2,r2) == 0) cycle ! ���� ����� �� ��������� � ���������, ������������ ���
            if (flag_status(c2,r2) == 2 .or. flag_status(c2,r2) == 3) cycle ! ���� ����� ��� ��� ��������� �����, ���� ������������ ���
            
            ! ��������� ���������� ��� ������
            waylength(c2,r2) = waylength(c1,r1) + sqrt( ( (c2-c1)*StepX )**2 + ( (r2-r1)*StepY )**2 )
            ! �������� ���� �������� ��� ������
            way_prev(c2,r2) = two_to_one(c1,r1)
            ! �������� ������ ��� ���������� ���������
            flag_status(c2,r2) = 3
            ! ��������� �������������� ������ � ������
            Select case (depr(c2,r2)) ! ��������� ������� � ����������� �� ����
            case(-2) ! ����� � ����� ������
                ! ��������� ����������� ����� ������ ����� �������� �� ��������� connect_depressions
                ! ���� TRUE, ������� ��������� �������������� ��� ������ �����
                ! ���� FALSE, ������ ��������� ����� �������������� ��������
                Select case (connect_depressions)
                case (.true.)
                    depr(c2,r2) = 1
                    flag_status(c2,r2) = 3
                    call add_cell_1(c2,r2,q1,q2)
                    
                    !!!call q%enqueue(Z(c2,r2), c2, r2)! ��������� � �������
                    
                case (.false.)
                    flag_status(c2,r2) = 2
                    if (Z1 > input_model(c2,r2) .and. input_model(c2,r2) > input_model(col_out,row_out)) then
                    Z1 = input_model(c2,r2)
                    c3 = c2; r3 = r2
                    endif
                end select
            case (-1) ! ����� � �������
                ! �������� ���������� ������ � ������
                flag_status(c2,r2) = 2
                q3 = q3 + 1
                list_borders(q3) = two_to_one(c2,r2)
                ! ���������� ������ ����������� ������ ����� ��������� �����. �� ����� ��������� ������, ��� �� �����
                if (Z1 > input_model(c2,r2) .and. input_model(c2,r2) > input_model(col_out,row_out)) then
                    Z1 = input_model(c2,r2)
                    c3 = c2; r3 = r2
                endif
            
            case (0) ! ����� � ������ �����
                continue
            case default ! ����� � ������ �����
                flag_status(c2,r2) = 3
                call add_cell_1(c2,r2,q1,q2)  
            end select
        
        enddo
        q2 = q2+1
    enddo
       
    ! ���� ������ ��� ������������ (Z0 � Z1) 
    !! �������� ��������� ������ ����� ������
    !!call decreasing_outlet()
    Z0 = output_model(col_out,row_out)
    ! �������� ������������ ��������� ������.
    !!call increasing_border(Z1)
    
    ! � ���������� ������� �����: �������������� ������������� ��������� ������; ���� ����������, ���� "����������" (������ ������ � ����������)
    ! ������ ����� ��������� (Z1) � �������� (Z0) ������ ������������
    
    ! ������ ����� � �������� ����������� (�� ������ � ������). ������������� ����� ��� ������ ������ (flag_status == 2) ������
    c1 = 0;r1 = 0; c2 = 0; r2 = 0
    q2 = 0
            
    do q2 = 1, q3
        call one_to_two(c1,r1,list_borders(q2))
            
        if (flag_status(c1,r1) /= 2) cycle ! ���� ����� �� ��������� � �������, ������������ �
        ! ������ � �����, �� ������� ���� ��������������� �����.            
        dZ = input_model(c1,r1) - Z0
        !dZ = Z1 - Z0
        length = waylength(c1,r1)
        do
            s = way_prev(c1,r1)
            call one_to_two(c2,r2,s)
            if (c2 == col_out .and. r2 == row_out) exit ! �������, ���� ����������� ����� �������� ������ ������
            Z_temp = Z0 + (dZ/length)*waylength(c2,r2)
            if (output_model(c2,r2) == nodata .or. output_model(c2,r2) > Z_temp) output_model(c2,r2) = Z_temp
            flag_status(c2,r2) = 1
            c1 = c2; r1 = r2
        enddo      
    enddo
    
    flag_status(col_out,row_out) = 1
    
    ! ���������� �����, ����� ������� �� �������� ���������� ����
    ! ���������� ������������ ���� �������������. ������������: ������ ����� ������ (Z0) � ������ ���������_���_������������ ����� (Z2)
    do i = 1, q1
        call one_to_two(c3,r3,list_cells(i))
        if (flag_status(c3,r3) /= 3) cycle
        dist1 = waylength(c3,r3) ! ������ ���� �� �������������� ������ �� ������
        c1 = c3; r1 = r3
        ! ������ ���� ���������� ��� ������������� ������ ���� �� �������   
        do
            s = way_prev(c1,r1)
            call one_to_two(c2,r2,s)
            if (flag_status(c2,r2) == -1) then ! ����, �������� �������, �� ����� �� ������-������, ��...
                Z2 = output_model(c2,r2) + 0.1
                dist2 = dist1
            else if (flag_status(c2,r2) == 1) then ! ����, �������� �������, �� ����� �� ����� � ��� ��������� �������, ��...
                Z2 = output_model(c2,r2) ! ������������� ����� �������� �������� ������ (������ ������ ����������� ����� �����)
                dist2 = waylength(c2,r2) ! � ���������� ���������� �� ������������� ����� �� ������
                exit ! � ������� �� �����
            else
                c1 = c2; r1 = r2 ! ���� ���, �� ������� ������
            endif
        enddo
        
        Z3 = Z0 + (Z2-Z0)/dist2*dist1       
        output_model(c3,r3) = Z3
        
        continue
        
    enddo

    ! ����������, ��. ������� output_model ��������� ������� � ������� ���������

    if(allocated(waylength)) deallocate(waylength)
    if(allocated(flag_status)) deallocate(flag_status)
    if(allocated(way_prev)) deallocate(way_prev)
    if(allocated(list_cells)) deallocate(list_cells)
    if(allocated(list_Z)) deallocate(list_Z)
    if(allocated(list_borders)) deallocate(list_borders)

contains

! ������������ ���������� ������

subroutine add_cell_1(c,r,q1,q2)

integer :: c,r,q1,q2 ! ������� � ������ ������, ����������� � ������, � ����� ����� �������� �������� � ������
real(4) :: Z_index
integer :: index
integer :: q
integer :: new_element ! ������, ��� ������� ����� ������� ����� �������. �� ��������� � � ����� ������
integer :: first, last
integer :: ncase
    
    index = two_to_one(c,r)

    q1 = q1+1
    new_element = q1 ! ����������, ��� ���, "�� ���������"
    
    if (q1 > 1) then ! ����������� �������� ���������� ������ �����, ����� � ������ ������ ������ ��������
        if (input_model(c,r) < list_Z(q1-1)) then 
            ! ���� ����� ������� ������ ����������, ��������� ��������� ������.
            ! (���� ����� ������� ������ ����������, ����� ������ ��� � �����)
            
            first = q2+1; last = q1-1
            ! ���� ����� �������� � ������. ���������� ����� (������, � �������)
            if (last /= first) then
                do while (last /= first)
                if (input_model(c,r) < list_Z( (last-first)/2 + first) ) then
                    last = (last-first)/2 + first
                else
                    first = (last-first)/2 + first + 1
                endif
                enddo
            endif
            new_element = first
        call move_cells_in_queue_1(new_element,q1)
        endif   
    endif        
    
    ! ���������� ����� ������� � ��� ������ � �� �������, ������� ���� ���������� (q, ��� �� new_element)
    list_cells(new_element) = index
    list_Z(new_element) = input_model(c,r)
end subroutine add_cell_1

subroutine move_cells_in_queue_1(first_cell,last_cell)
integer :: number
integer :: first_cell,last_cell

do number = last_cell, first_cell+1, -1
    list_cells(number) = list_cells(number-1)
    list_Z(number) = list_Z(number-1)
enddo

end subroutine move_cells_in_queue_1

! ������������ �������� ������ ����� ������
subroutine decreasing_outlet()
    real(4) :: Z0_temp
    real(4),parameter :: threshold = 0.0001
    
    Z0_temp = output_model(col_out,row_out)
    k = 0
    do k = 1,8
        c2 = col_out+kx(k); r2 = row_out+ky(k)
        
        if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
        if (depr(c2,r2) > 0) cycle ! �� ������������� ������ �� ���������
        if (output_model(c2,r2) + threshold < Z0_temp) then
            Z0_temp = output_model(c2,r2) + threshold
        endif
    enddo
    if (output_model(col_out,row_out) > Z0_temp) then
        output_model(col_out,row_out) = Z0_temp 
        print *, "Outlet point height decreased from ", input_model(col_out,row_out), " to ", Z0_temp
    else
        continue
        print *, "Outlet point height did not decreased"
    endif

end subroutine decreasing_outlet

! ������������ ���������� ��������� �����
subroutine increasing_border(Z1)
    real, allocatable :: matrix_Z1_temp(:,:)
    real(4),parameter :: threshold = 1.0
    real(4) :: Z1
    real(4) :: Z1_temp
    integer :: c,r
    
    ! ������������ �������� � ��� �����:
    ! 1) ��� ������ ����� �������, �� ������� ����� ���� ������������ (flag_status = 2), ������ ������������ ������������ �������� ������
    ! ���� ����� �������� ������ ������ (�� ��� ����� ���������!), �� � �������� �� ������������� (������������ �������� ����� ������������ ��������)
    ! 2) ����� ���� ������������ ����� (��� �����������!) ������ �����������;
    ! 3) ����������� ������ ����������� �� ������ ������������, ��� ������, ������� ���� � (�� ���� �������� ������ ����� ������!), �������������� �� � ��������    
    
    allocate(matrix_Z1_temp(Nx,Ny)) ! ������� ������������ �����
    matrix_Z1_temp(1:Nx,1:Ny) = nodata
    Z1_temp = Zmax ! �������������� �����, ���������� �������� ������
    c = 0; r = 0
    ! 1. ���������� ������������� ������ ����� (��� ������ ������)
    do c1 = 1,Nx
        do r1 = 1,Ny
            if (flag_status(c1,r1) /= 2) cycle
            matrix_Z1_temp(c1,r1) = input_model(c1,r1)
            if (depr(c1,r1) /= -1) cycle ! ��� ���� �����, ������� ��� ������ ��������, ��������� ������������ ������
            ! ��� �����, ������� ������� ��������, �������� ����������. ������������� ��� ���� ���� ���� �������
            do k = 1,8
                c2 = c1+kx(k); r2 = r1+ky(k)
                if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
                if (depr(c2,r2) > 0) cycle ! ������� �� ��������� ����������
                if (input_model(c2,r2) < input_model(c1,r1) .or. input_model(c2,r2) == nodata) cycle ! � ������ ����� ������ ������� ����������
                if (matrix_Z1_temp(c1,r1) < (input_model(c2,r2) - threshold) ) then
                    matrix_Z1_temp(c1,r1) = input_model(c2,r2) - threshold
                endif
            enddo
            ! 2. ���� �������� ����� ����������
            if (Z1_temp > matrix_Z1_temp(c1,r1)) then
                Z1_temp = matrix_Z1_temp(c1,r1); c = c1; r = r1
            endif
        enddo
    enddo
    
    ! ��������, ��������� �� ������
    if (Z1_temp == input_model(c,r)) then
        c = 0; r = 0
    endif
    
    Select case (c+r)
    case(0) ! ������ �� ���������
        print *, "Border point height did not increased"
    case default ! ������ ���������
        Z1 = Z1_temp
        do c1 = 1,Nx
            do r1 = 1,Ny
                if (flag_status(c1,r1) /= 2) cycle
                if (depr(c1,r1) /= -1) cycle ! ��, ��� �� ���� ������� (�� ����� ������!), ���������� �� �����
                if (output_model(c1,r1) > Z1) cycle ! ���� ������ ��� ���� Z1, ������ � �����������
                output_model(c1,r1) = Z1 ! ����������� ��� ����� ������ ������ �� ���������.
            enddo
        enddo
        print *, "Border point height increased from ", input_model(c,r), " to ", Z1
    end select
    
    deallocate(matrix_Z1_temp)
    
end subroutine increasing_border

integer function two_to_one(column,row)
integer :: row,column
two_to_one = (column-1)*Ny+row
end function two_to_one

subroutine one_to_two(column,row,index)
integer :: column,row,index
row = mod(index,Ny)
if (row == 0) then ! ���� �� ������ ����� ��������, ������ ����� ������������ �����������
    row = Ny
    index = index - Ny
endif
column = index/Ny + 1
end subroutine one_to_two

end subroutine filling_ASTAR