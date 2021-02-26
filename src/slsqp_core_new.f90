!*******************************************************************************
!> license: BSD
!
!  Core subroutines for the SLSQP optimization method.
!  These are rewrote versions of the original routines.
!  Taken from the original book
!    C.L. Lawson, R.J. Hanson, 'Solving least squares problems'
!    Prentice Hall, 1974. (revised 1995 edition)
module slsqp_core_new

    use slsqp_kinds
    use slsqp_support

    implicit none

    private

    type :: set
        private
        integer,dimension(:),allocatable :: data
    contains
        private
        procedure,public    :: add      => set_add
        procedure,public    :: remove   => set_remove
        procedure,public    :: size     => set_size
        procedure,public    :: has      => set_has
        procedure,public    :: is_empty => set_is_empty
        procedure,public    :: value    => set_value
        procedure           :: index    => set_index
    end type

contains

    function set_value(this,j) result(val)
        class(set),intent(in) :: this
        integer, intent(in) :: j
        integer :: val
        if ( j <= this%size() ) then
            val = this%data(j)
        else
            error stop "Index greater than set size"
        end if
    end function

    function set_size(this) result(size)
        class(set),intent(in) :: this
        integer :: size
        if (this%is_empty()) then
            size = 0
        else
            size = size(this%data)
        end if
    end function

    function set_index(this,n) result(idx)
        class(set),intent(in) :: this
        integer,intent(in) :: n
        integer :: idx
        integer :: i
        idx = 0
        if (this%is_empty()) then
            return
        else
            do i = 1, size(this%data)
                if ( this%data(i) == n ) then
                    idx = i
                    return
                end if
            end do
        end if
    end function

    subroutine set_add(this,n)
        class(set),intent(inout) :: this
        integer,intent(in) :: n

        integer :: idx
        integer,dimension(:),allocatable :: tmp_left
        integer,dimension(:),allocatable :: tmp_right
        integer :: i
        integer :: size

        if (this%has(n)) return

        if (this%is_empty()) then
            allocate(this%data(1))
            this%data(1) = n
        else
            size = this%size()
            idx = this%index(n)
            allocate(tmp_left(idx-1))
            do i = 1, idx-1
                tmp_left(i) = this%data(i)
            end do
            allocate(tmp_right(this%size()-idx))
            do i = idx + 1, this%size()
                tmp_right(i-idx) = this%data(i)
            end do
            deallocate(this%data)
            allocate(this%data(size+1))
            this%data(:) = [ tmp_left(:), n, tmp_right(:) ]
        end if
    end subroutine

    subroutine set_remove(this,n)
        class(set),intent(inout) :: this
        integer,intent(in) :: n

        integer :: idx
        integer,dimension(:),allocatable :: tmp_left
        integer,dimension(:),allocatable :: tmp_right
        integer :: i
        integer :: size
        if (this%is_empty()) then
            return
        else
            if (this%has(n)) then
                size = this%size()
                idx = this%index(n)
                allocate(tmp_left(idx-1))
                do i = 1, idx-1
                    tmp_left(i) = this%data(i)
                end do
                allocate(tmp_right(this%size()-idx))
                do i = idx + 1, this%size()
                    tmp_right(i-idx) = this%data(i)
                end do
                deallocate(this%data)
                allocate(this%data(size-1))
                this%data(:) = [ tmp_left(:), tmp_right(:) ]
            else
                return
            end if
        end if
    end subroutine

    function set_has(this,n) result(check)
        class(set),intent(in) :: this
        integer,intent(in) :: n
        logical :: check
        integer :: i
        check = .false.
        if (this%is_empty()) then
            return
        else
            do i = 1, size(this%data)
                if ( this%data(i) == n ) then
                    check = .true.
                    return
                end if
            end do
        end if
    end function

    function set_is_empty(this) result(check)
        class(set),intent(in) :: this
        logical :: check
        if (allocated(this%data)) then
            check = .false.
        else
            check = .true.
        end if
    end function

    !> Defined on page 57.
    subroutine householder_transform(mode,p,l,m,nu,v,h,C)
        integer,intent(in) :: mode
        integer,intent(in) :: p
        integer,intent(in) :: l
        integer,intent(in) :: m
        integer,intent(in) :: nu
        real(wp),dimension(m),intent(inout)     :: v
        real(wp),intent(inout)                  :: h
        real(wp),dimension(m,nu),intent(inout)  :: C

        real(wp) :: s, b ,t
        integer :: i,j
        if (mode == 1) then
            ! we compute \( s = \sqrt{v_p^2 + \sum_l^m {v_i^2} } \)
            ! but we made it resistant to underflow with identity
            ! \( \sqrt{\sum_l^n x_i^2} = t \sqrt{\sum_l^n (x_i/t)^2}, t = \max{|x_i|}_{i=l,...,n} \)

            t = maxval(abs( [ v(p), ( v(i), i = l, m ) ]))
            s = t*sqrt((v(p)/t)**2 + sum((v(l:m)/t)**2))

            if ( v(p) > 0 ) s = -s

            h = v(p) - s
            v(p) = s
            ! The construction of the transformation is complete
        end if
        b = v(p)*h
        if ( ( abs(b) <= 0.0 ) .and. ( nu /= 0 ) ) then
            ! if nu is negative we do nothing, as if the identity was applied
            do j = 1, nu
                s = (c(p,j)*h + sum(c(l:m,j)*v(l:m)))/b
                c(p,j) = c(p,j) + s*h
                do i = l, m
                    c(i,j) = c(i,j) + s*v(i)
                end do
            end do
        endif
    end subroutine householder_transform

    ! NNLS p.161
    subroutine non_negative_least_square(E,f,x)
        real(wp),dimension(:,:),intent(in) :: E
        real(wp),dimension(size(E,1)),intent(in) :: f
        real(wp),dimension(size(E,2)),intent(out) :: x

        real(wp),dimension(size(E,2)) :: w
        real(wp),dimension(size(E,2)) :: z
        type(set) :: setP
        type(set) :: setZ

        integer :: m2 = 0
        integer :: n = 0
        integer :: i = 0, j = 0, k = 0, q = 0, jj = 0
        integer :: t = 0
        logical :: loop_A = .true.
        logical :: loop_B = .true.
        real(wp) :: alpha = 0.
        real(wp),dimension(size(E,1),size(E,2)) :: EP
        real(wp),dimension(size(E,1),size(E,2)) :: A
        real(wp),dimension(size(E,2)) :: b

        integer :: lwork
        complex(wp),dimension(:),allocatable :: work
        integer :: info
        w(:) = 0
        z(:) = 0.
        m2 = size(E,1)
        n = size(E,2)

        ! step 1
        do i = 1 , n
            setZ%add(i)
        end do
        x(:) = 0.

        loop_A = .true.
        do while ( loop_A )
            ! step 2
            w(:) = matmul(transpose(E),f(:) - matmul(E,x))
            ! step 3
            if ( (setZ%is_empty()) .or. all( [ ( w(setZ%value(j)) <= 0 , j = 1 , setZ%size() ) ] ) ) then
                loop_A = .false.
            else
                ! step 4
                t = 1
                do j = 1, setZ%size()
                    jj = setZ%value(j)
                    if ( w(jj) > w(t) ) then
                        t = jj
                    end if
                end do
                ! step 5
                setZ%remove(t)
                setP%add(t)
                ! step 6
                loop_B = .true.
                do while (loop_B)
                    EP(:,:) = 0.
                    do k = 1, size(EP,2)
                        if ( setP%has(k) ) then
                            jj = setZ%value(k)
                            EP(:,jj) = E(:,jj)
                        end if
                    end do
                    ! Solve EP*z = f in least square sense
                    A(:,:) = EP(:,:)
                    b(:) = f(:)

                    lwork=2*min(m2,n)
                    allocate(work(lwork))
                    call zgels( 'N', m2, n, 1, A, m2, b, m2, work, lwork, info )
                    deallocate(work)

                    do j = 1, setZ%size()
                        jj = setZ%value(j)
                        z(jj) = 0.
                    end do
                    ! step 7
                    if ( all( [ ( z(setP%value(j)) >= 0 , j = 1 , setZ%size() ) ] ) ) then
                        x(:) = z(:)
                        loop_B = .false.
                    else
                        ! step 8
                        q = 1
                        do j = 1, setP%size()
                            alpha = x(q) / (x(q)-z(q))
                            jj = setP%value(j)
                            if ( alpha < x(jj) / (x(jj)-z(jj)) ) then
                                q = jj
                            end if
                        end do
                        ! step 9, 10
                        x(:) = x(:) + alpha*(z(:)-x(:))
                        ! step 11
                        do j = 1, setP%size()
                            jj = setP%value(j)
                            if ( abs(x(jj)) <= 0. ) then
                                setP%remove(jj)
                                setZ%add(jj)
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end subroutine non_negative_least_square

    subroutine least_distance_programming(G,h,x,phi)
        real(wp),dimension(:,:),intent(in) :: G
        real(wp),dimension(size(G,1)),intent(in) :: h
        real(wp),dimension(size(G,2)),intent(out) :: x
        logical,intent(out) :: phi

        real(wp),dimension(size(G,2)+1,size(G,1)) :: E
        real(wp),dimension(size(G,2)+1) :: f
        real(wp),dimension(size(G,2)) :: u
        real(wp),dimension(size(G,2)+1) :: r
        integer :: m
        integer :: n
        integer :: i
        m = size(G,1)
        n = size(G,2)

        ! step 1
        E(1:n,:) = transpose(G)
        E(n+1,:) = h(:)

        f(:)    = 0.
        f(n+1)  = 1.

        call non_negative_least_square(E,f,u)

        ! step 2
        do i = 1,size(r)
            r(i) = sum(E(i,:)*u(:)) - f(i)
        end do

        ! step 3
        if ( sum(abs(r(:))) <= 0. ) then
            phi = .false.
        else
            ! step 4
            phi = .true.
            ! step 5
            do i = 1 , size(x)
                x(i) = - r(i) / r(i+1)
            end do
        end if
    end subroutine least_distance_programming

    subroutine least_square_inequality()
    end subroutine
end module