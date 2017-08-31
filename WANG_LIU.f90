! Программа Wang_Liu — модифицированный вариант алгоритма Ванга-Лю для заполнения понижений
! Эта программа обрабатывает исходную ЦМР (Z), возвращая её «заполненный» вариант (Z_flat), а также матрицу с отмеченными понижениями и плоскими участками (depr)
    
    subroutine Wang_Liu(Z,Z_flat,depr,out_list,q_out,Zmax,Nx,Ny,NODATA)
    
    implicit none
    
    ! Переменные, передаваемые в подпрограмму
    real :: Z(Nx,Ny) ! Исходная ЦМР
    real :: Z_flat(Nx,Ny) ! Результат заполнения (под плоскость)
    integer :: depr(Nx,Ny) ! Идентификатор плоскостей на модели. 1 - плоскость, 0 - не плоскость.
    integer :: out_list(Nx*Ny) ! Список точек выхода (одномерный, сортировка по возрастанию высоты)
    integer :: q_out ! Итератор для списка точек выхода
    real :: Zmax ! Максимальная высота в пределах модели
    integer :: Nx,Ny ! Размерность сетки
    real :: NODATA ! Значение «нет данных»

    ! Внутренние переменные подпрограммы
    integer, allocatable :: mask(:,:) ! Массив, в котором будет отмечаться просмотренность/непросмотренность ячейки. 0 - не просмотрена, 1 - в очереди, 2 - просмотрена
    integer, allocatable :: list_cells(:) ! Общий список ячеек
    real, allocatable :: list_Z(:) ! Список высот ячеек. Параллелен списку list_cells. Необходим для удобного сравнения высоты (не надо каждый раз обращаться к исходной матрице)
    integer :: i,j ! Итераторы в циклах
    integer :: q1 ! Общее количество элементов в списке ячеек
    integer :: q2, q3 ! Итераторы в циклах (перебор ячеек в общем списке)
    integer :: q_max ! Общее количество точек для обработки. В начале работы программы равно количеству ячеек в матрице, затем снижается (за счёт значений "нет данных")
    integer :: index
    real :: Ztres
    integer :: percent_complete, percent_0
    integer :: c1,c2,r1,r2
    
    !! Скользящее окно 3х3 (поворот от направления "вправо" против часовой стрелки)
    !integer :: k
    !integer, dimension(8), parameter :: kx = [ 1, 1, 0,-1,-1,-1, 0, 1]
    !integer, dimension(8), parameter :: ky = [ 0, 1, 1, 1, 0,-1,-1,-1]
    
    ! Скользящее окно 3х3 (сначала просмотр прямых соседей, потом соседей по диагонали)
    integer :: k
    integer, dimension(8), parameter :: kx = [ 1, 0,-1, 0, 1,-1,-1, 1]
    integer, dimension(8), parameter :: ky = [ 0, 1, 0,-1, 1, 1,-1,-1]
    
! Начало работы
    ! Подготовка необходимых переменных и массивов.
    q_max = Nx*Ny
    Z_flat(1:Nx,1:Ny) = Z(1:Nx,1:Ny)
    
    depr(1:Nx,1:Ny) = 0
    allocate(mask(Nx,Ny))
    mask(1:Nx,1:Ny) = 0
    allocate(list_cells(Nx*Ny))
    allocate(list_Z(Nx*Ny))
    
    percent_0 = 0
    
    ! Определение первоначального списка точек
    q2 = 0
    call initialize_set()
    
    ! Обработка ЦМР
    q2 = 1
    q_out = 0
    do 
        if (q2 > q_max) exit
        call one_to_two(c1,r1,list_cells(q2)) ! Выбираем ячейку для просмотра на текущем шаге
    
        do k = 1,8 ! Идём в гости к соседям
            c2 = c1+kx(k); r2 = r1+ky(k)
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) cycle
            if (mask(c2,r2) /= 0) cycle ! Если у этого соседа мы уже были, переходим на следующую итерацию
            
            ! Добавляем новонайденного соседа в список
            call add_cell(c2,r2,q2)
                       
            ! Проверяем, есть ли в соседе понижение
            if (Z_flat(c2,r2) <= Z_flat(c1,r1)) then
            Z_flat(c2,r2) = Z_flat(c1,r1)
            depr(c2,r2) = 1
            
            if( depr(c1,r1) == 0) then ! Если точка, из которой мы пришли, не является понижением, маркируем её как точку выхода
                depr(c1,r1) = -2
                q_out = q_out + 1
                out_list(q_out) = two_to_one(c1,r1)
            endif
            
            endif
        enddo
        mask(c1,r1) = 2
        q2 = q2+1
        ! Отображение процента выполнения
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

! Подпрограмма формирования первоначального списка
subroutine initialize_set()
! Идея такая:
! Просматриваются все ячейки ЦМР. По умолчанию считается, что ячейку не следует добавлять в первоначальный список (decision = .false.)
! Однако, если ячейка находится на границе, или у неё сосед "нет данных" — решение меняется на положительное (decision = .true.)
! Ячейки со значениями "нет данных" в список не заносятся, а сразу убираются из рассмотрения
! В конце подпрограммы, если решение положительное, ячейка добавляется в первоначальный список

logical :: decision

q1 = 0
! Просматриваем по очереди все ячейки матрицы
do c1 = 1,Nx
    do r1 = 1,Ny
        decision = .false.
        ! Если в ячейке нет данных, исключаем её из рассмотрения насовсем
        if (Z(c1,r1) == nodata) then
            mask(c1,r1) = 2
            q_max = q_max - 1
            cycle
        endif
        ! Если в ячейке есть данные, присматриваемся к ней повнимательнее
        do k = 1,8
            c2 = c1+kx(k); r2 = r1+ky(k) ! Идем к соседу
            if (c2<1.or.c2>Nx.or.r2<1.or.r2>Ny) then ! Если сосед не существует вовсе, принимаем решение занести ячейку в список
                decision = .true.
                exit
            endif
            if (Z(c2,r2) == nodata) then ! Если сосед существует, но в нём "нет данных", принимаем решение занести ячейку в список
                decision = .true.
            endif
        enddo
        ! Добавление ячейки в список
        if (decision == .true.) then
            call add_cell(c1,r1,0)
        endif     
        
    enddo
enddo


end subroutine initialize_set

! Подпрограмма добавления элемента в список
subroutine add_cell(c,r,q2)

integer :: c,r,q2 ! Столбец и строка ячейки, добавляемой в список, а также номер текущего элемента в списке
real(4) :: Z_index
integer :: index
integer :: q
integer :: new_element ! Индекс, под которым будет записан новый элемент. По умолчанию — в конец списка
integer :: first, last
integer :: ncase
    
    index = two_to_one(c,r)

    q1 = q1+1
    new_element = q1 ! Собственно, вот оно, "по умолчанию"
    
    if (q1 > 1) then
        if (Z(c,r) < List_z(q1-1)) then ! Если новый элемент меньше последнего, запускаем процедуру сдвига.
            ! Если новый элемент больше последнего, сразу ставим его в конец
            first = q2+1; last = q1-1
            
            ! Попытка оптимизации
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
    
    ! Записываем новый элемент в оба списка в ту позицию, которая была определена (q, она же new_element)
    mask(c,r) = 1
    list_cells(new_element) = index
    list_Z(new_element) = Z(c,r)
end subroutine add_cell

! Подпрограмма перемещения элементов в списке
subroutine move_cells_in_queue(first_cell,last_cell)
integer :: number
integer :: first_cell,last_cell

do number = last_cell, first_cell+1, -1
    list_cells(number) = list_cells(number-1)
    list_Z(number) = list_Z(number-1)
enddo

end subroutine move_cells_in_queue

! Функция для конвертации двумерного индекса в одномерный
integer function two_to_one(column,row)
integer :: row,column
two_to_one = (column-1)*Ny+row
end function two_to_one

! Подпрограмма для восстановления двумерного индекса из одномерного
subroutine one_to_two(column,row,index)
integer :: column,row,index
row = mod(index,Ny)
if (row == 0) then ! Если не ввести такую поправку, индекс будет возвращаться неправильно
    row = Ny
    index = index - Ny
endif
column = index/Ny + 1    
end subroutine one_to_two

end subroutine Wang_Liu