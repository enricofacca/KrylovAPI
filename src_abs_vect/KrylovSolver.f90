module KrylovSolver  
  use AbstractOperators
  use DVectors
  implicit none
  private
  public :: pcg
  !! Derived type decscribing the operator
  !! y=(~A)^{-1} x
  !! with (~A)^{-1} is the approximate inverse via
  !! Preconditoned Conjugate Gradient
  type, public, extends(abstract_operator) :: pcg_solver
     !! set main controls 
     integer :: max_iterations
     double precision :: tolerance
     !! pointers for matrix and precontioner
     class(abstract_operator), pointer :: matrix
     class(abstract_operator), pointer :: preconditioner
     !! scratch vector
     class(dvector), allocatable :: aux(:)
   contains
     procedure, public, pass :: init => init_solver
     procedure, public, pass :: kill => kill_solver
     !! procedure overrided 
     procedure, public, pass :: apply => apply_pcg_solver
  end type pcg_solver
contains
  subroutine pcg(matrix,rhs,sol,ierr,&
       tolerance,max_iterations,preconditioner)
    use DVectors
    use AbstractOperators
    use IdentityOperators
    implicit none
    class(abstract_operator),           intent(inout) :: matrix
    class(dvector),                     intent(in   ) :: rhs
    class(dvector),                     intent(inout) :: sol
    integer,                            intent(inout) :: ierr
    double precision,                   intent(in   ) :: tolerance
    integer,                            intent(in   ) :: max_iterations
    class(abstract_operator), target,optional, intent(inout) :: preconditioner
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
    class(dvector), pointer :: aux(:)
    type(dvector), target, allocatable :: aux_local(:)
    class(dvector), pointer :: axp, pres,  resid, pk,scr
    double precision :: alpha, beta, presnorm,resnorm,bnorm
    double precision :: ptap
    double precision :: normres
    double precision :: tol,rort
    double precision :: inverse_residum_weight
    type(eye),target :: identity
    class(abstract_operator), pointer :: prec
 
   
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

    tol=tolerance
    max_iter=max_iterations
    
    

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
    allocate(aux_local(5))
    do i=1,5
       call aux_local(i)%init(nequ)
    end do
    aux=>aux_local

    axp          => aux(1)
    pres         => aux(2)
    resid        => aux(3)
    pk           => aux(4)
    scr          => aux(5)

    ! set tolerance
    tol = tolerance
    
    ! set max iterations number
    max_iter = max_iterations
    
    
    exit_test = .false.
    ! compute rhs norm and return zero solution
    bnorm = norm_dvector(rhs)

    if (bnorm<1.0d-16) then
       ! set solution and flag
       ierr=0
       do i=1,nequ
          call sol%set(i,0.0d0)
       end do

       ! free memory
       aux=>null()
       if (allocated(aux_local))then
          do i=1,6
             call aux_local(i)%kill()
          end do
          deallocate(aux_local)
       end if
       return
    end if
    
    

    ! calculate initial residual (res = rhs-M*sol)
    call matrix%apply(sol,resid,ierr)
    call axpby_dvector(1.0d0,rhs,-1.0d0,resid)
    resnorm = norm_dvector(resid)/bnorm
    
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
          beta = 0.0d0
       else
          beta = -dot_dvector(pres,axp)/ptap
       end if

       !  calculates p_{k+1}:=pres_{k}+beta_{k}*p_{k}
       call axpby_dvector(1.0d0,pres,beta,pk)
          
       !  calculates axp_{k+1}:= matrix * p_{k+1}
       call matrix%apply(pk,axp,ierr)
       
       !  calculates \alpha_k
       ptap  = dot_dvector(pk,axp)
       alpha = dot_dvector(pk,resid)/ptap
       
       !  calculates x_k+1 and r_k+1
       call axpby_dvector(alpha,pk,1.0d0,sol)
       call axpby_dvector(-alpha,axp,1.0d0,resid)

       !  compute residum
       resnorm = norm_dvector(resid)/bnorm
     
       ! checks
       exit_test = (&
            iter    .gt. max_iter .or.&
            resnorm .le. tol)

       if (iter.ge. max_iter) ierr = 1
    end do

    ! free memory
    aux=>null()
    if (allocated(aux_local))then
       do i=1,5
          call aux_local(i)%kill()
       end do
       deallocate(aux_local)
    end if

    
  end subroutine pcg
  
  subroutine init_solver(this,&
       matrix, max_iter, tolerance, precondtioner)
    class(pcg_solver),                intent(inout) :: this
    class(abstract_operator), target, intent(in   ) :: matrix
    integer,                          intent(in   ) :: max_iter
    double precision,                intent(in   ) :: tolerance
    class(abstract_operator), target, intent(in   ) :: precondtioner
    integer :: i

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
    allocate(this%aux(6))
    do i=1,6
       call this%aux(i)%init(this%nrow)
    end do
  end subroutine init_solver

  subroutine kill_solver(this)
    class(pcg_solver), intent(inout) :: this
    integer :: i
    this%max_iterations   = 0
    this% tolerance = 0.0d0
    !! the non-linearity will be approximately the tolerance of the
    !! pcg solver
    this%nonlinearity =  0.0d0

    this%matrix => null()
    this%preconditioner => null()

    !! free memory
    do i=1,6
       call this%aux(i)%kill()
    end do
    deallocate(this%aux)
  end subroutine kill_solver

  recursive subroutine apply_pcg_solver(this,vec_in,vec_out,ierr)
    implicit none
    class(pcg_solver), intent(inout) :: this
    class(dvector),    intent(in   ) :: vec_in
    class(dvector),    intent(inout) :: vec_out
    integer,           intent(inout) :: ierr

    call pcg(this%matrix, vec_in, vec_out, ierr,&
         tolerance=this%tolerance,&
         max_iterations=this%max_iterations,&
         preconditioner=this%preconditioner)

  end subroutine apply_pcg_solver
end module KrylovSolver
