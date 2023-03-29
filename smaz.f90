program grease
    implicit none
    integer :: NX, NY, s_max, i, j, approx_order, fff
    real(8) :: L, W, gamma, delta, HH, eps, dx, dy, F, U, findmoment, M
    real(8), allocatable :: p(:,:), x(:,:), y(:,:), h(:,:)
    
    open(1, file = 'input.txt')
    read(1,*) L 
    read(1,*) W 
    read(1,*) gamma 
    read(1,*) delta 
    read(1,*) U 
    read(1,*) HH 
    read(1,*) eps 
    read(1,*) NX 
    read(1,*) NY 
    read(1,*) s_max 
    read(1,*) approx_order 
    read(1,*) fff 
    close(1)
    
    open(2, file = 'optimalangle.txt')
    
    allocate(p(NX,NY), x(NX,NY), y(NX,NY), h(NX,NY))
    
    call mesh_making(NX, NY, x, y, L, W, dx, dy)
    
    if (fff == 1) then
        call Find_angle(NX, NY, dx, dy, x, y, h, p, L, HH, delta, s_max, approx_order, eps, U, gamma)
    end if
    
    call comp_h(NX, NY, dx, L, delta, gamma, HH, h)
    
    call bound_cond(NX, NY, p, approx_order)
    
    call solver(NX, NY, x, y, dx, dy, h, p, s_max, approx_order, eps, U)
    
    call output(NX, NY, x, y, h, p)
    
    call find_F(NX, NY, x, y, dx, dy, p, F)
    print *, 'F =', F
    
    write(2,*) '"Прогиб стрелы", "Угол наклона вкладыша"'
    write(2,*) delta, gamma
    close(2)
    
    deallocate(p, x, y)
    pause
    end program
    
    subroutine mesh_making(NX, NY, x, y, L, W, dx, dy) 
    implicit none
    integer :: i, j, NX, NY
    real(8) :: dx, dy, L, W
    real(8), dimension(NX, NY) :: x, y
    
    dx = L / (NX - 1)
    dy = W / 2.0 / (NY - 1)
    
    do i = 1, NX
        do j = 1, NY
        
            x(i,j) = (i - 1) * dx 
            y(i,j) = (j - 1) * dy 
            
        end do
    end do
    
    end subroutine
    
    subroutine comp_h(NX, NY, dx, L, delta, gamma, HH, h) 
    integer :: i, j, NX, NY
    real(8) :: dx, delta, gamma, L, HH
    real(8), dimension(NX,NY) :: h
        
    do i = 1, NX
        
        h(i,:) = HH + delta - delta * (1.0 - (4.0 / 3.1415926535**2.0) * ((i - 1) * dx - L / 2.0)**2.0) &
        - gamma * ((i - 1) * dx - L / 2.0)
        
    end do
    
    end subroutine
    
    subroutine bound_cond(NX, NY, p, approx_order)
    implicit none
    integer :: NX, NY, approx_order
    real(8), dimension(NX, NY) :: p
    
    p(1,:) = 0.0
    p(NX,:) = 0.0
    if (approx_order == 1) then
        p(:,1) = p(:,2)
    else if (approx_order == 2) then
        p(:,1) = 4.0 / 3.0 * p(:,2) - p(:,3) / 3.0
    end if
    p(:,NY) = 0.0
    
    end subroutine
    
    subroutine output(NX, NY, x, y, h, p) 
    implicit none
    integer :: i, j, NX, NY
    real(8), dimension(NX,NY) :: h, x, y, p, p_anlitic
    
    open(1, file = 'result.plt')
    write(1,*) 'variables = "x", "y", "h", "p", "p_anlitic"'
    write(1,*) 'zone i=', NY, ', j=', NX
    do i = 1, NX
        do j = 1, NY    
            write(1,*) x(i,j), y(i,j), h(i,j), p(i,j), p_anlitic(i,j)
        end do
    end do
    close(1)
    
    end subroutine
    
    subroutine find_F(NX, NY, x, y, dx, dy, p, F) 
    implicit none
    integer :: i, j, NX, NY
    real(8) :: F, dx, dy
    real(8), dimension(NX,NY) :: x, y, p
    F = 0.0
    do i = 1, NX-1
        do j = 1, NY-1
            F = F + dx * dy *(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1) - 4.0)/4.0    
  
        end do
    end do
    
    end subroutine
    
    subroutine solver(NX, NY, x, y, dx, dy, h, p, s_max, approx_order, eps, U) 
    implicit none
    integer :: s, i, j, s_max, NX, NY, approx_order
    real(8) :: U, dx, dy, eps, c1, c2, c3, c4, c5, c6, res, res1, F, mom, FindMoment
    real(8), dimension(NX, NY) :: h, x, y, p, p_old
    
    p(:,:) = 0.0
    !call find_F(NX, NY, x, y, dx, dy, p, F)
    !print *, 'F =', F
    !Mom = FindMoment(x, NX, NY, dx, dy, p)
    !print*, 'Moment =', Mom 
    !pause
    open(1, file = 'Res.plt')
    write(1,*) 'variables = "s", "Res", "F"'
    
    do s = 1, s_max
    
        p_old(:,:) = p(:,:)
        
        do i = 2, NX - 1
            do j = 2, NY - 1
                
                c1 = ((h(i+1,j) + h(i,j)) / 2.0)**3.0 / dx**2
                c2 = ((h(i-1,j) + h(i,j)) / 2.0)**3.0 / dx**2
                c3 = ((h(i,j+1) + h(i,j)) / 2.0)**3.0 / dy**2
                c4 = ((h(i,j-1) + h(i,j)) / 2.0)**3.0 / dy**2
                c5 = U * (h(i+1,j) - h(i-1,j)) / (2.0 * dx)
                c6 = c1 + c2 + c3 + c4
                
                p(i,j) = (c1 * p(i+1,j) + c2 * p(i-1,j) + c3 * p(i,j+1) + c4 * p(i,j-1) - c5) / c6
                
                call bound_cond(NX, NY, p, approx_order)
                
            end do
        end do
        
        call find_F(NX, NY, x, y, dx, dy, p, F) 
        
        res = maxval(abs(p - p_old))
        
        if (s == 1) then 
            res1 = res
        end if
        
        write(1,*) s, res / res1, F 
        if (mod(s, 100) == 0) then
            print *, 's =', s, 'Res =', res / res1
        end if
        
        if (res / res1 < eps) exit 
        
    end do
    
    close(1)
    
    end subroutine
    
    subroutine Find_angle(NX, NY, dx, dy, x, y, h, p, L, HH, delta, s_max, approx_order, eps, U, gamma) 
    integer :: NX, NY, o, OO, approx_order, s_max
    real(8) :: gamma1, gamma2, gamma, dx, dy, L, HH, delta, eps, U, M, M1, M2, FindMoment
    real(8), dimension(NX,NY) :: x, y, h, p
    
    OO = 100
    gamma1 = 0.0
    gamma2 = 6.0
    
    do o = 1, OO
        
        gamma = (gamma1 + gamma2) / 2.0
        
        if (abs(gamma2 - gamma1) < 1e-6) exit
        
        call comp_h(NX, NY, dx, L, delta, gamma, HH, h)
        call bound_cond(NX, NY, p, approx_order)
        call solver(NX, NY, x, y, dx, dy, h, p, s_max, approx_order, eps, U)
        M1 = FindMoment(x, NX, NY, dx, dy, p)
        
        call comp_h(NX, NY, dx, L, delta, gamma2, HH, h)
        call bound_cond(NX, NY, p, approx_order)
        call solver(NX, NY, x, y, dx, dy, h, p, s_max, approx_order, eps, U)
        M2 = FindMoment(x, NX, NY, dx, dy, p)
        
        if (M1 * M2 < 0.0) then
            gamma1 = gamma
        else
            gamma2 = gamma
        end if
        
    end do
    M = (M1+M2)/2
    print*, 'delta =', delta, 'angle =', gamma, 'M =', M
    
    end subroutine
    
 

    real(8) function FindMoment(x, NX, NY, dx, dy, p) 
    implicit none
    integer :: NX, NY, i, j
    real(8) :: dx, dy, L, R,  F
    real(8), dimension(NX,NY) :: x, p
    
    R=0.0
    L=0.0
    
    do i = 1, NX-1 
        F = 0.0
        do j = 1, NY-1   
            F = F + dx * dy * (p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1) - 4.0) / 4.0    
        enddo        
        if (i < 0.5*(NX+1)) then
            R = F*(0.5-(x(i+1,1)+x(i,1))/2.0) + R
        else
            L = F*((x(i,1)+x(i+1,1))/2.0-0.5) + L
        end if
    end do 
    FindMoment = R - L
    
    end function

    subroutine analytical_solution(NX, NY, p_anlitic, x, h, L, U, mu)
        implicit none
        integer :: i, j, NX, NY
        real(8) :: L, U, mu
        real(8), dimension(NX, NY) :: p_anlitic, x, h
        
        do i = 1, NX
            do j = 1, NY
                p_anlitic(i,j) = 6.0 * mu * U * L / ((h(1,j)**2 - h(NX,j)**2) * (h(1,j) - h(i,j)) * (h(i,j) - h(NX,j)) / h(i,j)**2)
            end do
        end do
        
        end subroutine
