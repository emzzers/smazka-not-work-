program cr_pr
    implicit none
    real(8) W, L, U, pi, dx, dy, delta, HH
    real(8) eps, F, eps1, gamma, a, b, Sum_M, fa, fc, c, mu
    integer i, j, NI, NJ, s_max, approx_order, fff
    real(8), allocatable :: x(:,:), y(:,:), P(:,:), h(:,:), p_analytic(:,:)

    open(1, file='input.txt')
    open(2, file='Res.plt')
    open(3, file='result.plt')
    pi = 4.0 * atan(1.0)
    read(1,*) L
    read(1,*) W
    read(1,*) gamma
    read(1,*) delta
    read(1,*) U
    read(1,*) HH
    read(1,*) eps
    read(1,*) NI
    read(1,*) NJ
    read(1,*) s_max
    read(1,*) approx_order
    read(1,*) fff
    close(1)
    allocate(x(NI,NJ), y(NI,NJ), P(NI,NJ), h(NI,NJ))
    a = 0.0
    b = 2.0 * (1.0 + delta)
    call xy_array(dx, dy, x, y, L, W, NI, NJ)
    call h_array(gamma, L, delta, pi, x, h, NI, NJ, HH)
    call p_init(NI, NJ, P, approx_order)
    call p_F_M(NI, NJ, x, dx, dy, h, P, s_max, F, eps, U, Sum_M)
    print *, F
    eps1 = 1.0e-6
    fa = F
    do while (abs(a-b) > eps1)
        c = (a+b)/2.0
        call h_array(c, L, delta, pi, x, h, NI, NJ, HH)
        call p_init(NI, NJ, P, approx_order)
        call p_F_M(NI, NJ, x, dx, dy, h, P, s_max, F, eps, U, Sum_M)
        fc = F
        if ((fa*fc) < 0) then
            b = c
        else
            a = c
            fa = fc
        endif
        print *, fc
        !pause
    enddo

    gamma = c
    Sum_M = fc
    print *, gamma, F, Sum_M
    p_analytic = 6.0 * mu * U * L / ((h(1,:)**2 - h(NI,:)**2) * (h(1,:) - h(:,1)) * (h(:,1) - h(NI,:)) / h(:,:)**2)

    call output(NI, NJ, X, Y, P, p_analytic, h)

    deallocate(x, y, P, h, p_analytic)
    close(1)
    close(2)
    close(3)
    pause
end program

subroutine xy_array(dx, dy, x, y, L, W, NI, NJ)
    implicit none
    integer NI, NJ, i, j
    real(8) dx, dy, L, W
    real(8), dimension(NI,NJ) :: x, y
    dx=L/(NI-1.0) 
    dy=W/2.0/(NJ-1.0) 
    x(1,1:NJ)=0.0
    y(1:NI,1)=0.0
    do j=1,NJ
    do i=2,NI
        x(i,j)=x(i-1,j)+dx 
    enddo
    enddo
    do j=2,NJ
    do i=1,NI
        y(i,j)=y(i,j-1)+dy 
    enddo
    enddo
    end subroutine

    subroutine h_array (gamma,L,delta,pi,x,h,NI,NJ,HH)
    implicit none
    integer i,NI,NJ
    real(8) gamma,L,delta,pi, HH
    real(8),dimension(NI,NJ):: x, h
    do i=1,NI
        h(i,:) = HH + delta - delta * (1.0 - (4.0 / pi**2.0) * (x(i,1) - L / 2.0)**2.0) &
    - gamma * (x(i,1) - L / 2.0)
    enddo
    end subroutine
    
    subroutine p_init(NI,NJ,P, approx_order)
        implicit none
        integer NI,NJ, approx_order
        real(8), dimension (NI,NJ):: P
        p(1,:) = 0.0
        p(NI,:) = 0.0
        if (approx_order == 1) then
            p(:,1) = p(:,2)
        else if (approx_order == 2) then
            p(:,1) = 4.0 / 3.0 * p(:,2) - p(:,3) / 3.0
        end if
        p(:,NJ) = 0.0
    end subroutine

    subroutine p_F_M(NI, NJ, x, dx, dy, h, p, s_max, F, eps, U, Sum_M)
    implicit none
    integer s, i, j, s_max, NI, NJ, approx_order
    real(8)  U, dx, dy, eps, c1, c2, c3, c4, c5, c6, res, res1, F, F0, M1, M2, Sum_M
    real(8), dimension (NI,NJ):: x, p, p_old, h
    p(:,:) = 0.0

    write(2,*) 'variables = "s", "Res", "F"'
    do s = 1, s_max
        p_old(:,:) = p(:,:)
     
        do i=2,NI-1
            do j=2,NJ-1
            c1 = ((h(i+1,j) + h(i,j)) / 2.0)**3.0 / dx**2
            c2 = ((h(i-1,j) + h(i,j)) / 2.0)**3.0 / dx**2
            c3 = ((h(i,j+1) + h(i,j)) / 2.0)**3.0 / dy**2
            c4 = ((h(i,j-1) + h(i,j)) / 2.0)**3.0 / dy**2
            c5 = U * (h(i+1,j) - h(i-1,j)) / (2.0 * dx)
            c6 = c1 + c2 + c3 + c4
            
            p(i,j) = (c1 * p(i+1,j) + c2 * p(i-1,j) + c3 * p(i,j+1) + c4 * p(i,j-1) - c5) / c6
            call p_init(NI,NJ,P, approx_order)
            enddo
    
        enddo
        F=0.0
        M1=0.0
        M2=0.0
        do i=1,NI-1
            F0=0.0
            do j=1,NJ-1
                F0=F0+dx*dy*(P(i,j)+P(i+1,j)+P(i,j+1)+P(i+1,j+1)-4.0)/4.0
                F=F+dx*dy*(P(i,j)+P(i+1,j)+P(i,j+1)+P(i+1,j+1)-4.0)/4.0
            enddo
        
            if (i < 0.5*(NI+1)) then
                M1=M1+(0.5-(x(i+1,1)+x(i,1))/2.0)*F0
            else
                M2=M2+((x(i,1)+x(i+1,1))/2.0-0.5)*F0
            endif
        enddo
        Sum_M=M1-M2

        res = maxval(abs(p - p_old))
    
        if (s == 1) then
            res1 = res
        end if
        
        write(2,*) s, res / res1, F 
        
        if (mod(s, 10) == 0) then
            print *, 's =', s, 'Res =', res / res1
        end if
        
        if (res / res1 < eps) exit 

    enddo
    
    
    end subroutine

    subroutine output(NI,NJ,X,Y,P,p_anlitic,h)
        implicit none
        integer NI,NJ
        real(8), dimension (NI,NJ):: X,Y,P,p_anlitic, h
      
        Write(3,*) 'VARIABLES = "X", "Y", "P", "p_anlitic", "h"'
        Write(3,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
        Write(3,'(100E25.16)') x(1:NI,1:NJ) 
        Write(3,'(100E25.16)') y(1:NI,1:NJ)    
        Write(3,'(100E25.16)') h(1:NI,1:NJ)     
        Write(3,'(100E25.16)') p(1:NI,1:NJ)
        Write(3,'(100E25.16)') p_anlitic(1:NI,1:NJ)    
        end subroutine
    