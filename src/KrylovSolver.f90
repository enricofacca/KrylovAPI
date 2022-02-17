module KrylovSolver  
  use AbstractOperators
  private
  public :: pcg
  !! Derived type decscribing the operator
  !! y=(~A)^{-1} x
  !! with (~A)^{-1} is the approximate inverse via
  !! Preconditoned Conjugate Gradient
  type, public, extends(abstract_operator) :: pcg_solver
     !! set main controls 
     integer :: max_iterations
     real(kind=double) :: tolerance
     !! pointers for matrix and precontioner
     class(abstract_operator), pointer :: matrix
     class(abstract_operator), pointer :: preconditioner
     !! scratch vector
     real(kind=double), allocatable :: aux(:)
   contains
     procedure, public, pass :: init => init_solver
     procedure, public, pass :: kill => kill_solver
     !! procedure overrided 
     procedure, public, pass :: apply => apply_pcg_solver
  end type pcg_solver
contains
  subroutine pcg(matrix,rhs,sol,ierr,&
       tolerance,max_iterations,preconditioner,aux_passed)
    use AbstractOperators
    use SimpleVectors
    implicit none
    class(abstract_operator),           intent(inout) :: matrix
    real(kind=double),                  intent(in   ) :: rhs(matrix%ncol)
    real(kind=double),                  intent(inout) :: sol(matrix%nrow)
    integer,                            intent(inout) :: ierr
    real(kind=double),        optional, intent(in   ) :: tolerance
    integer,                  optional, intent(in   ) :: max_iterations
    class(abstract_operator), target,optional, intent(inout) :: preconditioner
    real(kind=double), target,optional, intent(inout) :: aux_passed(6*matrix%ncol)
    !
    ! local vars
    !
    ! logical
    logical :: exit_test
    ! string
    character(len=256) :: msg
    ! integer
    integer :: i
    integer :: lun_err, lun_out 
    integer :: iprt,max_iter
    integer :: nequ,dim_ker,ndir=0
    integer :: iort
    integer :: info_prec
    integer :: iter
    ! real
    real(kind=double), pointer :: aux(:)
    real(kind=double), target, allocatable :: aux_local(:)
    real(kind=double), pointer :: axp(:), pres(:), ainv(:), resid(:), pk(:),scr(:)
    real(kind=double) :: alpha, beta, presnorm,resnorm,bnorm
    real(kind=double) :: ptap
    real(kind=double) :: normres
    real(kind=double) :: tol,rort
    real(kind=double) :: inverse_residum_weight
    !
    type(eye), target :: identity
    class(abstract_operator), pointer :: prec => null()

    !! checks
    if  (.not. matrix%is_symmetric)  then
       ierr = 998
       return
    end if        
    if  (matrix%nrow .ne. matrix%ncol )  then
       ierr = 999
       return
    end if           
    nequ=matrix%nrow

    if (present(tolerance)) then
       tol=tolerance
    else
       tol=1.0d-6
    end if

     if (present(max_iterations)) then
       max_iter=max_iterations
    else
       max_iter=200
    end if
    
    
    

    !! handle preconditing
    if (present(preconditioner)) then
       if  (.not. preconditioner%is_symmetric)  then
          ierr = 999
          return
       end if
       if  ( preconditioner%nrow .ne. preconditioner%ncol )  then
          ierr = 999
          return
       end if
       prec => preconditioner
    else
       call identity%init(nequ)
       prec => identity
    end if

    ! allocate memory if required and set pointers for scratch arrays
    if (present(aux_passed)) then
       aux=>aux_passed
    else
       allocate(aux_local(6*nequ))
       aux=>aux_local
    end if
    axp          => aux(       1 :   nequ)
    pres         => aux(  nequ+1 : 2*nequ)
    ainv         => aux(2*nequ+1 : 3*nequ)
    resid        => aux(3*nequ+1 : 4*nequ)
    pk           => aux(4*nequ+1 : 5*nequ)
    scr          => aux(5*nequ+1 : 6*nequ)

    ! set tolerance
    if (present (tolerance)) then
       tol = tolerance
    else
       tol = 1.0d-6
    end if

    ! set max iterations number
    if (present (max_iterations)) then
       max_iter = max_iterations
    else
       max_iter = 100
    end if

    
    exit_test = .false.
    ! compute rhs norm and break cycle 
    bnorm = d_norm(nequ,rhs)
    if (bnorm<1.0d-16) then
       ierr=0
       sol=zero
       ! free memory
       aux=>null()
       if (allocated(aux_local)) deallocate(aux_local)
       return
    end if
    
    

    ! calculate initial residual (res = rhs-M*sol)
    call matrix%apply(sol,resid,ierr)
    resid = rhs - resid
    resnorm = d_norm(nequ,resid)/bnorm
    
    !
    ! cycle
    !
    iter=0
    do while (.not. exit_test)
       iter = iter + 1   
       ! compute  pres = PREC  (r_{k+1})
       call prec%apply(resid,pres,ierr)

       !  calculates beta_k
       if (iter.eq.1) then
          beta = zero
       else
          beta = -d_scalar_product(nequ,pres,axp)/ptap
       end if

       !  calculates p_{k+1}:=pres_{k}+beta_{k}*p_{k}
       call d_axpby(nequ,one,pres,beta,pk)
          
       !  calculates axp_{k+1}:= matrix * p_{k+1}
       call matrix%apply(pk,axp,ierr)
       
       !  calculates \alpha_k
       ptap  = d_scalar_product(nequ,pk,axp)
       alpha = d_scalar_product(nequ,pk,resid)/ptap
       
       !  calculates x_k+1 and r_k+1
       call d_axpby(nequ,alpha,pk,one,sol)
       call d_axpby(nequ,-alpha,axp,one,resid)

       !  compute residum
       resnorm = d_norm(nequ,resid)/bnorm
     
       ! checks
       exit_test = (&
            iter    .gt. max_iter .or.&
            resnorm .le. tol)

       if (iter.ge. max_iter) ierr = 1
    end do

    ! free memory
    aux=>null()
    if (allocated(aux_local)) deallocate(aux_local)

    
  end subroutine pcg
  
  subroutine init_solver(this,&
       matrix, max_iter, tolerance, precondtioner)
    class(pcg_solver),                intent(inout) :: this
    class(abstract_operator), target, intent(in   ) :: matrix
    integer,                          intent(in   ) :: max_iter
    real(kind=double),                intent(in   ) :: tolerance
    class(abstract_operator), target, intent(in   ) :: precondtioner

    !! the set domain/codomain operator dimension
    this%nrow = matrix%nrow
    this%ncol = matrix%ncol

    !! the inverse of a symmetric matrix is symmetric but the
    !! approximiate inverse it is not in general.
    !! We set that that the operator is symmettric if
    !! if the require tolerance is below 1e-13
    if (tolerance < 1e-13) then
       this%is_symmetric=.True.
    else
       this%is_symmetric=.False.
    end if
       
    !! the non-linearity will be (approximately) the tolerance of the
    !! pcg solver
    this%nonlinearity =  tolerance


    !! set controls
    this%max_iterations  = max_iter
    this%tolerance = tolerance

    !! set pointer
    this%matrix => matrix
    this%preconditioner => precondtioner

    !! set work arrays
    allocate(this%aux(6*matrix%nrow))
    aux = zero
  end subroutine init_solver

  subroutine kill_solver(this)
    class(pcg_solver), intent(inout) :: this

    this%max_iterations   = 0
    this% tolerance = zero
    !! the non-linearity will be approximately the tolerance of the
    !! pcg solver
    this%nonlinearity =  zero

    this%matrix => null()
    this%preconditioner => null()

    !! free memory
    deallocate(this%aux)
  end subroutine kill_solver

  subroutine apply_pcg_solver(this,vec_in,vec_out,ierr)
    implicit none
    class(pcg_solver),  intent(inout) :: this
    real(kind=double),  intent(in   ) :: vec_in(this%ncol)
    real(kind=double),  intent(inout) :: vec_out(this%nrow)
    integer,                   intent(inout) :: ierr

    call pcg(this%matrix, vec_in, vec_out, ierr,&
         tolerance=this%tolerance,&
         max_iterations=this%max_iterations,&
         preconditioner=this%preconditioner,&
         aux_passed=this%aux)

  end subroutine apply_pcg_solver
end module KrylovSolver
