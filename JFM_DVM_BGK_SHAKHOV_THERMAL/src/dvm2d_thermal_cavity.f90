program dvm2d_thermal_cavity
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: eps = 1.0e-300_dp
  character(len=*), parameter :: solver_version = '2026-07-21-jfm-thermal-dvm-v3'

  integer :: nx, ny, nv, steps, min_steps, log_every, save_every
  integer :: projection_max_iter
  real(dp) :: vmin, vmax, kn, t_hot, t_cold, pr, cfl, tol, floor_val
  real(dp) :: shakhov_weight_min, shakhov_weight_max, projection_tol
  character(len=32) :: model
  character(len=256) :: out_prefix
  logical :: converged

  real(dp) :: dx, dy, dt, vmax_abs, time_now, final_res
  real(dp) :: initial_mass, final_mass, mass_drift_relative
  real(dp) :: max_target_conservation_error, max_physical_wall_mass_flux
  integer(kind=8) :: negative_target_values, clipped_updated_values
  integer(kind=8) :: shakhov_limited_target_values
  real(dp) :: shakhov_raw_weight_min, shakhov_raw_weight_max
  real(dp) :: shakhov_limited_h_base_mass, shakhov_total_h_base_mass
  real(dp) :: shakhov_limited_b_base_content, shakhov_total_b_base_content
  integer(kind=8) :: positive_projection_failures
  integer :: max_projection_iterations_used
  integer :: it, last_iter

  real(dp), allocatable :: cu(:), cv(:), wu(:), wv(:), w(:,:)
  real(dp), allocatable :: x(:,:), y(:,:)
  real(dp), allocatable :: h(:,:,:,:), b(:,:,:,:), hadv(:,:,:,:), badv(:,:,:,:)
  real(dp), allocatable :: rho(:,:), ux(:,:), uy(:,:), temp(:,:), press(:,:), qx(:,:), qy(:,:)
  real(dp), allocatable :: thetax(:,:), thetay(:,:), thetaz(:,:), sigxx(:,:), sigyy(:,:), sigxy(:,:)
  real(dp), allocatable :: m3x(:,:), m3y(:,:), m4x(:,:), m4y(:,:), sx(:,:), sy(:,:), kx(:,:), ky(:,:)
  real(dp), allocatable :: rho_old(:,:), ux_old(:,:), uy_old(:,:), temp_old(:,:)

  namelist /params/ nx, ny, nv, vmin, vmax, kn, t_hot, t_cold, pr, &
       cfl, tol, floor_val, shakhov_weight_min, shakhov_weight_max, &
       projection_max_iter, projection_tol, &
       steps, min_steps, log_every, save_every, model, out_prefix

  call set_defaults()
  call read_input()
  call validate_input()

  if (mod(nv-1,4) /= 0) then
     nv = (nv/4)*4 + 1
     write(*,*) 'Adjusted nv to Boole-compatible value:', nv
  end if

  allocate(cu(nv), cv(nv), wu(nv), wv(nv), w(nv,nv))
  call boole_grid(nv, vmin, vmax, cu, wu)
  call boole_grid(nv, vmin, vmax, cv, wv)
  call make_weight(nv, wu, wv, w)

  ! One quarter of the paper domain is solved: x,y in [0,1/2].  The planes
  ! x=0 and y=0 are exact specular symmetry boundaries.  The physical right
  ! and top walls are diffuse at T_cold and T_hot, respectively.
  dx = 0.5_dp / real(nx, dp)
  dy = 0.5_dp / real(ny, dp)
  vmax_abs = max(maxval(abs(cu)), maxval(abs(cv)))
  dt = cfl / (vmax_abs/dx + vmax_abs/dy + eps)

  allocate(x(nx,ny), y(nx,ny))
  allocate(h(nx,ny,nv,nv), b(nx,ny,nv,nv), hadv(nx,ny,nv,nv), badv(nx,ny,nv,nv))
  allocate(rho(nx,ny), ux(nx,ny), uy(nx,ny), temp(nx,ny), press(nx,ny), qx(nx,ny), qy(nx,ny))
  allocate(thetax(nx,ny), thetay(nx,ny), thetaz(nx,ny), sigxx(nx,ny), sigyy(nx,ny), sigxy(nx,ny))
  allocate(m3x(nx,ny), m3y(nx,ny), m4x(nx,ny), m4y(nx,ny), sx(nx,ny), sy(nx,ny), kx(nx,ny), ky(nx,ny))
  allocate(rho_old(nx,ny), ux_old(nx,ny), uy_old(nx,ny), temp_old(nx,ny))

  call init_geometry()
  call init_distribution()
  call compute_moments_all(h, b)
  initial_mass = 4.0_dp*sum(rho)*dx*dy
  max_target_conservation_error = 0.0_dp
  max_physical_wall_mass_flux = 0.0_dp
  negative_target_values = 0_8
  clipped_updated_values = 0_8
  shakhov_limited_target_values = 0_8
  shakhov_raw_weight_min = huge(1.0_dp)
  shakhov_raw_weight_max = -huge(1.0_dp)
  shakhov_limited_h_base_mass = 0.0_dp
  shakhov_total_h_base_mass = 0.0_dp
  shakhov_limited_b_base_content = 0.0_dp
  shakhov_total_b_base_content = 0.0_dp
  positive_projection_failures = 0_8
  max_projection_iterations_used = 0
  rho_old = rho; ux_old = ux; uy_old = uy; temp_old = temp

  open(unit=20, file=trim(out_prefix)//'.hst', status='replace', action='write')
  write(20,*) 'VARIABLES = iter, sim_time, dt, res_rho, res_u, res_v, res_T'

  write(*,*) 'JFM stationary thermal-cavity DVM run'
  write(*,*) 'solver_version=', solver_version
  write(*,*) 'quarter nx,ny,nv=', nx, ny, nv
  write(*,*) 'model=', trim(model), ' paper Kn=', kn, ' RT=', t_cold/t_hot, ' dt=', dt
  write(*,*) 'collision time: tau=2*Kn/(rho*sqrt(pi*T))'

  converged = .false.
  final_res = huge(1.0_dp)
  last_iter = 0

  do it = 1, steps
     last_iter = it
     call advect_step()
     call collide_step()
     call compute_moments_all(h, b)
     time_now = real(it, dp)*dt

     if (mod(it, log_every) == 0 .or. it == 1) then
        call residual_and_log(it, time_now, final_res)
        if (it >= min_steps .and. final_res < tol) then
           converged = .true.
           write(*,*) 'CONVERGED at iter=', it, ' final_res=', final_res
           exit
        end if
     end if

     if (save_every > 0 .and. mod(it, save_every) == 0) then
        call write_all_outputs(trim(out_prefix)//'_checkpoint_'//trim(int_to_str(it)))
     end if
  end do

  close(20)
  call compute_boundary_flux_diagnostic(max_physical_wall_mass_flux)
  final_mass = 4.0_dp*sum(rho)*dx*dy
  mass_drift_relative = abs(final_mass-initial_mass)/(abs(initial_mass)+eps)
  call write_all_outputs(trim(out_prefix))
  call write_metadata(last_iter, real(last_iter,dp)*dt, converged, final_res)
  write(*,*) 'DONE. final iter=', last_iter, ' time=', real(last_iter,dp)*dt, ' final_res=', final_res

contains

  subroutine set_defaults()
    nx = 51; ny = 51; nv = 61
    vmin = -5.0_dp; vmax = 5.0_dp
    kn = 30.0_dp; t_hot = 1.0_dp; t_cold = 0.2_dp
    pr = 2.0_dp/3.0_dp
    cfl = 0.45_dp; tol = 1.0e-8_dp; floor_val = 0.0_dp
    shakhov_weight_min = 0.0_dp; shakhov_weight_max = 2.0_dp
    projection_max_iter = 12; projection_tol = 1.0e-12_dp
    steps = 120000; min_steps = 20000; log_every = 200; save_every = 0
    model = 'shakhov'
    out_prefix = 'cavity_dvm'
  end subroutine set_defaults

  subroutine read_input()
    integer :: ios
    open(unit=10, file='dvm_cavity.in', status='old', action='read', iostat=ios)
    if (ios == 0) then
       read(10, nml=params)
       close(10)
       write(*,*) 'Read input file: dvm_cavity.in'
    else
       write(*,*) 'No dvm_cavity.in found. Using defaults.'
    end if
  end subroutine read_input

  subroutine validate_input()
    character(len=32) :: model_clean
    model_clean = trim(adjustl(model))
    if (model_clean /= 'bgk' .and. model_clean /= 'shakhov') then
       write(*,*) 'ERROR: model must be bgk or shakhov; got ', trim(model)
       error stop 2
    end if
    model = model_clean
    if (nx < 8 .or. ny < 8 .or. nv < 9) error stop 'nx, ny, or nv is too small'
    if (vmin >= 0.0_dp .or. vmax <= 0.0_dp .or. abs(vmax+vmin) > 1.0e-12_dp) &
         error stop 'velocity interval must be symmetric about zero'
    if (kn <= 0.0_dp) error stop 'paper Kn must be positive'
    if (t_hot <= 0.0_dp .or. t_cold <= 0.0_dp .or. t_cold >= t_hot) &
         error stop 'require 0 < t_cold < t_hot'
    if (abs(t_hot-1.0_dp) > 1.0e-12_dp) &
         error stop 'paper scaling requires t_hot=1'
    if (pr <= 0.0_dp .or. pr > 1.0_dp) error stop 'Pr must be in (0,1]'
    if (shakhov_weight_min < 0.0_dp .or. &
         shakhov_weight_max <= shakhov_weight_min) &
         error stop 'invalid Shakhov weight limiter'
    if (projection_max_iter < 2 .or. projection_tol <= 0.0_dp) &
         error stop 'invalid positive projection controls'
    if (steps < 1 .or. min_steps < 0 .or. log_every < 1) error stop 'invalid iteration controls'
  end subroutine validate_input

  subroutine boole_grid(n, a, bnd, xg, wg)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, bnd
    real(dp), intent(out) :: xg(n), wg(n)
    real(dp) :: hh
    integer :: i, i0
    real(dp), dimension(5) :: coeff
    hh = (bnd-a)/real(n-1,dp)
    do i=1,n
       xg(i) = a + real(i-1,dp)*hh
       wg(i) = 0.0_dp
    end do
    coeff = (/7.0_dp, 32.0_dp, 12.0_dp, 32.0_dp, 7.0_dp/) * (2.0_dp*hh/45.0_dp)
    do i0 = 1, n-1, 4
       wg(i0:i0+4) = wg(i0:i0+4) + coeff
    end do
  end subroutine boole_grid

  subroutine make_weight(n, wu, wv, ww)
    integer, intent(in) :: n
    real(dp), intent(in) :: wu(n), wv(n)
    real(dp), intent(out) :: ww(n,n)
    integer :: i,j
    do j=1,n
       do i=1,n
          ww(i,j) = wu(i)*wv(j)
       end do
    end do
  end subroutine make_weight

  subroutine init_geometry()
    integer :: i,j
    do j=1,ny
       do i=1,nx
          x(i,j) = (real(i,dp)-0.5_dp)*dx
          y(i,j) = (real(j,dp)-0.5_dp)*dy
       end do
    end do
  end subroutine init_geometry

  subroutine init_distribution()
    integer :: i,j
    do j=1,ny
       do i=1,nx
          call set_maxwellian_cell(1.0_dp, 0.0_dp, 0.0_dp, &
               0.5_dp*(t_hot+t_cold), h(i,j,:,:), b(i,j,:,:))
       end do
    end do
  end subroutine init_distribution

  subroutine set_maxwellian_cell(r0, u0, v0, T0, hh, bb)
    real(dp), intent(in) :: r0, u0, v0, T0
    real(dp), intent(out) :: hh(nv,nv), bb(nv,nv)
    integer :: ku, kv
    real(dp) :: mass, val, scale
    mass = 0.0_dp
    do kv=1,nv
       do ku=1,nv
          val = r0/(pi*T0) * exp(-((cu(ku)-u0)**2 + (cv(kv)-v0)**2)/T0)
          hh(ku,kv) = val
          mass = mass + val*w(ku,kv)
       end do
    end do
    scale = r0/(mass + eps)
    do kv=1,nv
       do ku=1,nv
          hh(ku,kv) = scale*hh(ku,kv)
          bb(ku,kv) = 0.5_dp*T0*hh(ku,kv)
       end do
    end do
  end subroutine set_maxwellian_cell

  real(dp) function wall_unit_M(ku, kv, uw, vw, Tw)
    integer, intent(in) :: ku, kv
    real(dp), intent(in) :: uw, vw, Tw
    wall_unit_M = 1.0_dp/(pi*Tw) * exp(-((cu(ku)-uw)**2 + (cv(kv)-vw)**2)/Tw)
  end function wall_unit_M

  real(dp) function diffuse_rho_wall(hadj, side, uw, vw, Tw)
    real(dp), intent(in) :: hadj(nv,nv)
    character(len=*), intent(in) :: side
    real(dp), intent(in) :: uw, vw, Tw
    integer :: ku, kv
    real(dp) :: numer, denom, M
    numer = 0.0_dp; denom = 0.0_dp
    select case(trim(side))
    case('left')
       do kv=1,nv; do ku=1,nv
          if (cu(ku) > 0.0_dp) then
             M = wall_unit_M(ku,kv,uw,vw,Tw)
             denom = denom + cu(ku)*M*w(ku,kv)
          else if (cu(ku) < 0.0_dp) then
             numer = numer - cu(ku)*hadj(ku,kv)*w(ku,kv)
          end if
       end do; end do
    case('right')
       do kv=1,nv; do ku=1,nv
          if (cu(ku) < 0.0_dp) then
             M = wall_unit_M(ku,kv,uw,vw,Tw)
             denom = denom + cu(ku)*M*w(ku,kv)
          else if (cu(ku) > 0.0_dp) then
             numer = numer - cu(ku)*hadj(ku,kv)*w(ku,kv)
          end if
       end do; end do
    case('bottom')
       do kv=1,nv; do ku=1,nv
          if (cv(kv) > 0.0_dp) then
             M = wall_unit_M(ku,kv,uw,vw,Tw)
             denom = denom + cv(kv)*M*w(ku,kv)
          else if (cv(kv) < 0.0_dp) then
             numer = numer - cv(kv)*hadj(ku,kv)*w(ku,kv)
          end if
       end do; end do
    case('top')
       do kv=1,nv; do ku=1,nv
          if (cv(kv) < 0.0_dp) then
             M = wall_unit_M(ku,kv,uw,vw,Tw)
             denom = denom + cv(kv)*M*w(ku,kv)
          else if (cv(kv) > 0.0_dp) then
             numer = numer - cv(kv)*hadj(ku,kv)*w(ku,kv)
          end if
       end do; end do
    end select
    diffuse_rho_wall = max(numer/(denom + eps), 1.0e-12_dp)
  end function diffuse_rho_wall

  subroutine advect_step()
    real(dp), allocatable :: rhoR(:), rhoT(:)
    integer :: i,j,ku,kv,kur,kvr
    real(dp) :: fL, fR, fB, fT, gL, gR, gB, gT
    real(dp) :: MR, MT
    allocate(rhoR(ny), rhoT(nx))

    do j=1,ny
       rhoR(j) = diffuse_rho_wall(h(nx,j,:,:), 'right', 0.0_dp, 0.0_dp, t_cold)
    end do
    do i=1,nx
       rhoT(i) = diffuse_rho_wall(h(i,ny,:,:), 'top', 0.0_dp, 0.0_dp, t_hot)
    end do

    !$omp parallel do collapse(4) private(i,j,ku,kv,kur,kvr,fL,fR,fB,fT,gL,gR,gB,gT,MR,MT)
    do kv=1,nv
       do ku=1,nv
          do j=1,ny
             do i=1,nx
                kur = nv + 1 - ku
                kvr = nv + 1 - kv

                ! x=0 symmetry plane: exact specular reflection.
                if (i == 1) then
                   if (cu(ku) >= 0.0_dp) then
                      fL = cu(ku)*h(i,j,kur,kv)
                      gL = cu(ku)*b(i,j,kur,kv)
                   else
                      fL = cu(ku)*h(i,j,ku,kv)
                      gL = cu(ku)*b(i,j,ku,kv)
                   end if
                else
                   if (cu(ku) >= 0.0_dp) then
                      fL = cu(ku)*h(i-1,j,ku,kv)
                      gL = cu(ku)*b(i-1,j,ku,kv)
                   else
                      fL = cu(ku)*h(i,j,ku,kv)
                      gL = cu(ku)*b(i,j,ku,kv)
                   end if
                end if

                ! x=1/2 physical cold diffuse wall.
                if (i == nx) then
                   if (cu(ku) >= 0.0_dp) then
                      fR = cu(ku)*h(i,j,ku,kv)
                      gR = cu(ku)*b(i,j,ku,kv)
                   else
                      MR = rhoR(j)*wall_unit_M(ku,kv,0.0_dp,0.0_dp,t_cold)
                      fR = cu(ku)*MR
                      gR = cu(ku)*(0.5_dp*t_cold*MR)
                   end if
                else
                   if (cu(ku) >= 0.0_dp) then
                      fR = cu(ku)*h(i,j,ku,kv)
                      gR = cu(ku)*b(i,j,ku,kv)
                   else
                      fR = cu(ku)*h(i+1,j,ku,kv)
                      gR = cu(ku)*b(i+1,j,ku,kv)
                   end if
                end if

                ! y=0 symmetry plane: exact specular reflection.
                if (j == 1) then
                   if (cv(kv) >= 0.0_dp) then
                      fB = cv(kv)*h(i,j,ku,kvr)
                      gB = cv(kv)*b(i,j,ku,kvr)
                   else
                      fB = cv(kv)*h(i,j,ku,kv)
                      gB = cv(kv)*b(i,j,ku,kv)
                   end if
                else
                   if (cv(kv) >= 0.0_dp) then
                      fB = cv(kv)*h(i,j-1,ku,kv)
                      gB = cv(kv)*b(i,j-1,ku,kv)
                   else
                      fB = cv(kv)*h(i,j,ku,kv)
                      gB = cv(kv)*b(i,j,ku,kv)
                   end if
                end if

                ! y=1/2 physical hot diffuse wall.
                if (j == ny) then
                   if (cv(kv) >= 0.0_dp) then
                      fT = cv(kv)*h(i,j,ku,kv)
                      gT = cv(kv)*b(i,j,ku,kv)
                   else
                      MT = rhoT(i)*wall_unit_M(ku,kv,0.0_dp,0.0_dp,t_hot)
                      fT = cv(kv)*MT
                      gT = cv(kv)*(0.5_dp*t_hot*MT)
                   end if
                else
                   if (cv(kv) >= 0.0_dp) then
                      fT = cv(kv)*h(i,j,ku,kv)
                      gT = cv(kv)*b(i,j,ku,kv)
                   else
                      fT = cv(kv)*h(i,j+1,ku,kv)
                      gT = cv(kv)*b(i,j+1,ku,kv)
                   end if
                end if

                hadv(i,j,ku,kv) = max(h(i,j,ku,kv) - dt/dx*(fR-fL) - dt/dy*(fT-fB), floor_val)
                badv(i,j,ku,kv) = max(b(i,j,ku,kv) - dt/dx*(gR-gL) - dt/dy*(gT-gB), floor_val)
             end do
          end do
       end do
    end do
    !$omp end parallel do

    deallocate(rhoR, rhoT)
  end subroutine advect_step

  subroutine collide_step()
    integer :: i,j,ku,kv
    integer :: kk,ll,proj_it,used_proj_iters
    real(dp) :: tau, ratio, Hbase, Bbase, Hraw, Braw, Htarget, Btarget
    real(dp) :: corr_h, corr_b, raw_corr_h, raw_corr_b
    real(dp) :: step_raw_min, step_raw_max, exponent_h, exponent_b
    real(dp) :: cx, cy, c2, theta, dotq, Acoef, pcell, updated_h, updated_b
    real(dp) :: psi_h(4), psi_b(4), amat(4,4), rhs(4), delta(4)
    real(dp) :: lambda(4), mtarg(4), desired(4), cres(4)
    real(dp) :: scale_cons, cell_cons_error, step_maxerr, projection_error
    real(dp) :: max_increment, step_length
    real(dp) :: step_limited_h_mass, step_total_h_mass
    real(dp) :: step_limited_b_content, step_total_b_content
    integer(kind=8) :: step_negative, step_clipped, step_limited, step_projection_failures
    integer :: step_max_projection_iterations
    logical :: solve_ok, projection_ok

    call compute_moments_all(hadv, badv)

    step_maxerr = 0.0_dp
    step_negative = 0_8
    step_clipped = 0_8
    step_limited = 0_8
    step_projection_failures = 0_8
    step_max_projection_iterations = 0
    step_limited_h_mass = 0.0_dp
    step_total_h_mass = 0.0_dp
    step_limited_b_content = 0.0_dp
    step_total_b_content = 0.0_dp
    step_raw_min = huge(1.0_dp)
    step_raw_max = -huge(1.0_dp)
    !$omp parallel do collapse(2) &
    !$omp& private(i,j,ku,kv,kk,ll,proj_it,used_proj_iters,tau,ratio,Hbase,Bbase,Hraw,Braw, &
    !$omp& Htarget,Btarget,corr_h,corr_b,raw_corr_h,raw_corr_b,exponent_h,exponent_b, &
    !$omp& cx,cy,c2,theta,dotq,Acoef,pcell,updated_h,updated_b,psi_h,psi_b,amat,rhs,delta, &
    !$omp& lambda,mtarg,desired,cres,scale_cons,cell_cons_error,projection_error, &
    !$omp& max_increment,step_length,solve_ok,projection_ok) &
    !$omp& reduction(max:step_maxerr,step_raw_max,step_max_projection_iterations) &
    !$omp& reduction(min:step_raw_min) &
    !$omp& reduction(+:step_negative,step_clipped,step_limited,step_projection_failures, &
    !$omp& step_limited_h_mass,step_total_h_mass,step_limited_b_content,step_total_b_content)
    do j=1,ny
       do i=1,nx
          ! This is exactly the model equation in the paper:
          !   streaming = (1/Kn) * [rho*sqrt(pi*T)/2] * (target-f),
          ! hence tau = 2*Kn/[rho*sqrt(pi*T)].
          tau = 2.0_dp*kn/((rho(i,j)+eps)*sqrt(pi*max(temp(i,j),1.0e-10_dp)))
          tau = max(tau, 1.0e-12_dp)
          ratio = dt/tau
          theta = 0.5_dp*max(temp(i,j),1.0e-10_dp)
          pcell = 0.5_dp*rho(i,j)*max(temp(i,j),1.0e-10_dp)
          desired(1) = rho(i,j)
          desired(2) = rho(i,j)*ux(i,j)
          desired(3) = rho(i,j)*uy(i,j)
          desired(4) = rho(i,j)*(0.75_dp*temp(i,j) + &
               0.5_dp*(ux(i,j)**2+uy(i,j)**2))

          ! Positive conservative projection.  The limited BGK/Shakhov target
          ! is exponentially tilted by four Lagrange multipliers.  Newton's
          ! method matches discrete mass, two momenta and total energy while
          ! preserving non-negativity at every velocity node.
          lambda = 0.0_dp
          projection_ok = .false.
          do proj_it=1,projection_max_iter
             amat = 0.0_dp
             mtarg = 0.0_dp
             do kv=1,nv
                do ku=1,nv
                   cx = cu(ku)-ux(i,j)
                   cy = cv(kv)-uy(i,j)
                   c2 = cx*cx + cy*cy
                   Hbase = rho(i,j)/(pi*temp(i,j)) * exp(-c2/temp(i,j))
                   Bbase = 0.5_dp*temp(i,j)*Hbase
                   Hraw = Hbase
                   Braw = Bbase
                   if (trim(model) == 'shakhov') then
                      dotq = cx*qx(i,j) + cy*qy(i,j)
                      Acoef = (1.0_dp-pr)*dotq/(5.0_dp*pcell*theta + eps)
                      raw_corr_h = 1.0_dp + Acoef*(c2/theta - 4.0_dp)
                      raw_corr_b = 1.0_dp + Acoef*(c2/theta - 2.0_dp)
                      if (proj_it == 1) then
                         step_total_h_mass = step_total_h_mass + Hbase*w(ku,kv)
                         step_total_b_content = step_total_b_content + Bbase*w(ku,kv)
                         step_raw_min = min(step_raw_min,raw_corr_h,raw_corr_b)
                         step_raw_max = max(step_raw_max,raw_corr_h,raw_corr_b)
                         if (raw_corr_h < shakhov_weight_min .or. &
                              raw_corr_h > shakhov_weight_max) then
                            step_limited = step_limited + 1_8
                            step_limited_h_mass = step_limited_h_mass + Hbase*w(ku,kv)
                         end if
                         if (raw_corr_b < shakhov_weight_min .or. &
                              raw_corr_b > shakhov_weight_max) then
                            step_limited = step_limited + 1_8
                            step_limited_b_content = step_limited_b_content + Bbase*w(ku,kv)
                         end if
                      end if
                      corr_h = min(max(raw_corr_h,shakhov_weight_min),shakhov_weight_max)
                      corr_b = min(max(raw_corr_b,shakhov_weight_min),shakhov_weight_max)
                      Hraw = Hbase*corr_h
                      Braw = Bbase*corr_b
                   end if

                   psi_h = (/1.0_dp, cu(ku), cv(kv), 0.5_dp*(cu(ku)**2+cv(kv)**2)/)
                   psi_b = (/0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp/)
                   exponent_h = max(-80.0_dp,min(80.0_dp,dot_product(lambda,psi_h)))
                   exponent_b = max(-80.0_dp,min(80.0_dp,0.5_dp*lambda(4)))
                   Htarget = Hraw*exp(exponent_h)
                   Btarget = Braw*exp(exponent_b)

                   mtarg(1) = mtarg(1) + Htarget*w(ku,kv)
                   mtarg(2) = mtarg(2) + cu(ku)*Htarget*w(ku,kv)
                   mtarg(3) = mtarg(3) + cv(kv)*Htarget*w(ku,kv)
                   mtarg(4) = mtarg(4) + 0.5_dp*((cu(ku)**2+cv(kv)**2)*Htarget+Btarget)*w(ku,kv)
                   do kk=1,4
                      do ll=1,4
                         amat(kk,ll) = amat(kk,ll) + w(ku,kv)*( &
                              Htarget*psi_h(kk)*psi_h(ll) + Btarget*psi_b(kk)*psi_b(ll))
                      end do
                   end do
                end do
             end do

             rhs = desired-mtarg
             scale_cons = max(rho(i,j),abs(desired(4)),1.0e-14_dp)
             projection_error = maxval(abs(rhs))/scale_cons
             if (projection_error <= projection_tol) then
                projection_ok = .true.
                exit
             end if
             call solve_4x4(amat, rhs, delta, solve_ok)
             if (.not. solve_ok) exit
             max_increment = abs(delta(1)) + vmax_abs*(abs(delta(2))+abs(delta(3))) + &
                  vmax_abs*vmax_abs*abs(delta(4))
             max_increment = max(max_increment,0.5_dp*abs(delta(4)))
             step_length = min(1.0_dp,2.0_dp/(max_increment+eps))
             lambda = lambda + step_length*delta
          end do
          used_proj_iters = min(proj_it,projection_max_iter)

          ! Rebuild the final positive target, verify its conserved moments,
          ! and take the implicit relaxation update.
          mtarg = 0.0_dp
          do kv=1,nv
             do ku=1,nv
                cx = cu(ku)-ux(i,j)
                cy = cv(kv)-uy(i,j)
                c2 = cx*cx + cy*cy
                Hbase = rho(i,j)/(pi*temp(i,j)) * exp(-c2/temp(i,j))
                Bbase = 0.5_dp*temp(i,j)*Hbase
                Hraw = Hbase
                Braw = Bbase
                if (trim(model) == 'shakhov') then
                   dotq = cx*qx(i,j) + cy*qy(i,j)
                   Acoef = (1.0_dp-pr)*dotq/(5.0_dp*pcell*theta + eps)
                   raw_corr_h = 1.0_dp + Acoef*(c2/theta - 4.0_dp)
                   raw_corr_b = 1.0_dp + Acoef*(c2/theta - 2.0_dp)
                   corr_h = min(max(raw_corr_h,shakhov_weight_min),shakhov_weight_max)
                   corr_b = min(max(raw_corr_b,shakhov_weight_min),shakhov_weight_max)
                   Hraw = Hbase*corr_h
                   Braw = Bbase*corr_b
                end if
                psi_h = (/1.0_dp, cu(ku), cv(kv), 0.5_dp*(cu(ku)**2+cv(kv)**2)/)
                exponent_h = max(-80.0_dp,min(80.0_dp,dot_product(lambda,psi_h)))
                exponent_b = max(-80.0_dp,min(80.0_dp,0.5_dp*lambda(4)))
                Htarget = Hraw*exp(exponent_h)
                Btarget = Braw*exp(exponent_b)
                mtarg(1) = mtarg(1) + Htarget*w(ku,kv)
                mtarg(2) = mtarg(2) + cu(ku)*Htarget*w(ku,kv)
                mtarg(3) = mtarg(3) + cv(kv)*Htarget*w(ku,kv)
                mtarg(4) = mtarg(4) + 0.5_dp*((cu(ku)**2+cv(kv)**2)*Htarget+Btarget)*w(ku,kv)
                if (Htarget < 0.0_dp) step_negative = step_negative + 1_8
                if (Btarget < 0.0_dp) step_negative = step_negative + 1_8
                updated_h = (hadv(i,j,ku,kv) + ratio*Htarget)/(1.0_dp+ratio)
                updated_b = (badv(i,j,ku,kv) + ratio*Btarget)/(1.0_dp+ratio)
                if (updated_h < floor_val) step_clipped = step_clipped + 1_8
                if (updated_b < floor_val) step_clipped = step_clipped + 1_8
                h(i,j,ku,kv) = max(updated_h, floor_val)
                b(i,j,ku,kv) = max(updated_b, floor_val)
             end do
          end do
          cres = mtarg-desired
          scale_cons = max(rho(i,j),abs(desired(4)),1.0e-14_dp)
          cell_cons_error = maxval(abs(cres))/scale_cons
          if (cell_cons_error <= projection_tol) projection_ok = .true.
          if (.not. projection_ok) step_projection_failures = step_projection_failures + 1_8
          step_maxerr = max(step_maxerr,cell_cons_error)
          step_max_projection_iterations = max(step_max_projection_iterations,used_proj_iters)
       end do
    end do
    !$omp end parallel do
    max_target_conservation_error = max(max_target_conservation_error,step_maxerr)
    negative_target_values = negative_target_values + step_negative
    clipped_updated_values = clipped_updated_values + step_clipped
    if (trim(model) == 'shakhov') then
       shakhov_limited_target_values = shakhov_limited_target_values + step_limited
       shakhov_raw_weight_min = min(shakhov_raw_weight_min,step_raw_min)
       shakhov_raw_weight_max = max(shakhov_raw_weight_max,step_raw_max)
       shakhov_limited_h_base_mass = shakhov_limited_h_base_mass + step_limited_h_mass
       shakhov_total_h_base_mass = shakhov_total_h_base_mass + step_total_h_mass
       shakhov_limited_b_base_content = shakhov_limited_b_base_content + step_limited_b_content
       shakhov_total_b_base_content = shakhov_total_b_base_content + step_total_b_content
    end if
    positive_projection_failures = positive_projection_failures + step_projection_failures
    max_projection_iterations_used = max(max_projection_iterations_used,step_max_projection_iterations)
  end subroutine collide_step

  subroutine solve_4x4(ain, bin, solution, ok)
    real(dp), intent(in) :: ain(4,4), bin(4)
    real(dp), intent(out) :: solution(4)
    logical, intent(out) :: ok
    real(dp) :: aa(4,4), bb(4), rowtmp(4), tmp, factor, pivabs
    integer :: i,j,k,piv
    aa = ain
    bb = bin
    ok = .true.
    do k=1,4
       piv = k
       pivabs = abs(aa(k,k))
       do i=k+1,4
          if (abs(aa(i,k)) > pivabs) then
             piv = i
             pivabs = abs(aa(i,k))
          end if
       end do
       if (pivabs < 1.0e-30_dp) then
          ok = .false.
          solution = 0.0_dp
          return
       end if
       if (piv /= k) then
          rowtmp = aa(k,:)
          aa(k,:) = aa(piv,:)
          aa(piv,:) = rowtmp
          tmp = bb(k); bb(k) = bb(piv); bb(piv) = tmp
       end if
       do i=k+1,4
          factor = aa(i,k)/aa(k,k)
          aa(i,k) = 0.0_dp
          do j=k+1,4
             aa(i,j) = aa(i,j)-factor*aa(k,j)
          end do
          bb(i) = bb(i)-factor*bb(k)
       end do
    end do
    solution = 0.0_dp
    do i=4,1,-1
       tmp = bb(i)
       do j=i+1,4
          tmp = tmp-aa(i,j)*solution(j)
       end do
       solution(i) = tmp/aa(i,i)
    end do
  end subroutine solve_4x4

  subroutine compute_moments_all(hh, bb)
    real(dp), intent(in) :: hh(nx,ny,nv,nv), bb(nx,ny,nv,nv)
    integer :: i,j,ku,kv
    real(dp) :: mass, mx, my, ener, cx, cy, c2int
    !$omp parallel do collapse(2) private(i,j,ku,kv,mass,mx,my,ener,cx,cy,c2int)
    do j=1,ny
       do i=1,nx
          mass = 0.0_dp; mx = 0.0_dp; my = 0.0_dp; ener = 0.0_dp
          do kv=1,nv
             do ku=1,nv
                mass = mass + hh(i,j,ku,kv)*w(ku,kv)
                mx = mx + cu(ku)*hh(i,j,ku,kv)*w(ku,kv)
                my = my + cv(kv)*hh(i,j,ku,kv)*w(ku,kv)
                ener = ener + 0.5_dp*((cu(ku)**2+cv(kv)**2)*hh(i,j,ku,kv)+bb(i,j,ku,kv))*w(ku,kv)
             end do
          end do
          rho(i,j) = max(mass, 1.0e-14_dp)
          ux(i,j) = mx/(rho(i,j)+eps)
          uy(i,j) = my/(rho(i,j)+eps)
          temp(i,j) = (4.0_dp/3.0_dp)*(ener/(rho(i,j)+eps) - 0.5_dp*(ux(i,j)**2+uy(i,j)**2))
          temp(i,j) = max(temp(i,j), 1.0e-10_dp)
          ! Paper pressure is scaled by rho0*R*T_hot.  Since molecular
          ! velocity is scaled by sqrt(2*R*T_hot), this is twice the raw
          ! velocity-space normal momentum flux.
          press(i,j) = rho(i,j)*temp(i,j)
          qx(i,j)=0.0_dp; qy(i,j)=0.0_dp
          thetax(i,j)=0.0_dp; thetay(i,j)=0.0_dp; thetaz(i,j)=0.0_dp
          sigxx(i,j)=0.0_dp; sigyy(i,j)=0.0_dp; sigxy(i,j)=0.0_dp
          m3x(i,j)=0.0_dp; m3y(i,j)=0.0_dp; m4x(i,j)=0.0_dp; m4y(i,j)=0.0_dp
          do kv=1,nv
             do ku=1,nv
                cx = cu(ku)-ux(i,j)
                cy = cv(kv)-uy(i,j)
                thetax(i,j) = thetax(i,j) + cx*cx*hh(i,j,ku,kv)*w(ku,kv)
                thetay(i,j) = thetay(i,j) + cy*cy*hh(i,j,ku,kv)*w(ku,kv)
                thetaz(i,j) = thetaz(i,j) + bb(i,j,ku,kv)*w(ku,kv)
                sigxy(i,j) = sigxy(i,j) + cx*cy*hh(i,j,ku,kv)*w(ku,kv)
                c2int = (cx*cx+cy*cy)*hh(i,j,ku,kv) + bb(i,j,ku,kv)
                qx(i,j) = qx(i,j) + 0.5_dp*cx*c2int*w(ku,kv)
                qy(i,j) = qy(i,j) + 0.5_dp*cy*c2int*w(ku,kv)
                m3x(i,j) = m3x(i,j) + cx**3*hh(i,j,ku,kv)*w(ku,kv)
                m3y(i,j) = m3y(i,j) + cy**3*hh(i,j,ku,kv)*w(ku,kv)
                m4x(i,j) = m4x(i,j) + cx**4*hh(i,j,ku,kv)*w(ku,kv)
                m4y(i,j) = m4y(i,j) + cy**4*hh(i,j,ku,kv)*w(ku,kv)
             end do
          end do
          thetax(i,j) = thetax(i,j)/(rho(i,j)+eps)
          thetay(i,j) = thetay(i,j)/(rho(i,j)+eps)
          thetaz(i,j) = thetaz(i,j)/(rho(i,j)+eps)
          sigxx(i,j) = 2.0_dp*rho(i,j)*thetax(i,j) - press(i,j)
          sigyy(i,j) = 2.0_dp*rho(i,j)*thetay(i,j) - press(i,j)
          sigxy(i,j) = 2.0_dp*sigxy(i,j)
          sx(i,j) = m3x(i,j)/(rho(i,j)*max(thetax(i,j),eps)**1.5_dp + eps)
          sy(i,j) = m3y(i,j)/(rho(i,j)*max(thetay(i,j),eps)**1.5_dp + eps)
          kx(i,j) = m4x(i,j)/(rho(i,j)*max(thetax(i,j),eps)**2.0_dp + eps)
          ky(i,j) = m4y(i,j)/(rho(i,j)*max(thetay(i,j),eps)**2.0_dp + eps)
       end do
    end do
    !$omp end parallel do
  end subroutine compute_moments_all

  subroutine residual_and_log(iter, tnow, resmax)
    integer, intent(in) :: iter
    real(dp), intent(in) :: tnow
    real(dp), intent(out) :: resmax
    real(dp) :: rrho, ru, rv, rT
    rrho = maxval(abs(rho-rho_old))/(maxval(abs(rho))+eps)
    ru = maxval(abs(ux-ux_old))/(maxval(abs(ux))+eps)
    rv = maxval(abs(uy-uy_old))/(maxval(abs(uy))+eps)
    rT = maxval(abs(temp-temp_old))/(maxval(abs(temp))+eps)
    resmax = max(max(rrho,ru), max(rv,rT))
    write(*,'(A,I9,A,ES13.5,A,ES13.5,A,4ES13.5)') 'iter=', iter, ' time=', tnow, ' dt=', dt, ' res=', rrho, ru, rv, rT
    write(20,'(I12,6ES18.8)') iter, tnow, dt, rrho, ru, rv, rT
    rho_old = rho; ux_old = ux; uy_old = uy; temp_old = temp
  end subroutine residual_and_log

  subroutine write_block(unit, A)
    integer, intent(in) :: unit
    real(dp), intent(in) :: A(nx,ny)
    integer :: i,j,count
    count = 0
    do j=1,ny
       do i=1,nx
          write(unit,'(ES23.16,2X)', advance='no') A(i,j)
          count = count + 1
          if (mod(count,6) == 0) write(unit,*)
       end do
    end do
    if (mod(count,6) /= 0) write(unit,*)
  end subroutine write_block

  subroutine write_vblock(unit, A)
    integer, intent(in) :: unit
    real(dp), intent(in) :: A(nv,nv)
    integer :: ku,kv,count
    count = 0
    do kv=1,nv
       do ku=1,nv
          write(unit,'(ES23.16,2X)', advance='no') A(ku,kv)
          count = count + 1
          if (mod(count,6) == 0) write(unit,*)
       end do
    end do
    if (mod(count,6) /= 0) write(unit,*)
  end subroutine write_vblock

  subroutine write_all_outputs(prefix)
    character(len=*), intent(in) :: prefix
    integer :: unit, is(6), js(6)
    character(len=24) :: names(6)
    real(dp) :: speed(nx,ny)
    call compute_moments_all(h,b)
    open(newunit=unit, file=trim(prefix)//'_quarter.dat', status='replace', action='write')
    write(unit,*) 'VARIABLES = X, Y, RHO, U, V, T, P, QX_PAPER, QY_PAPER'
    write(unit,*) 'ZONE I = ', nx, ', J = ', ny, ', DATAPACKING=BLOCK'
    call write_block(unit,x); call write_block(unit,y); call write_block(unit,rho); call write_block(unit,ux)
    call write_block(unit,uy); call write_block(unit,temp); call write_block(unit,press)
    call write_block(unit,2.0_dp*qx); call write_block(unit,2.0_dp*qy)
    close(unit)

    open(newunit=unit, file=trim(prefix)//'_quarter_moments.dat', status='replace', action='write')
    write(unit,*) 'VARIABLES = X, Y, RHO, U, V, T, P, QX_PAPER, QY_PAPER, THETAX, THETAY, THETAZ, SIGXX, SIGYY, SIGXY, M3X, M3Y, SX, SY, M4X, M4Y, KX, KY'
    write(unit,*) 'ZONE I = ', nx, ', J = ', ny, ', DATAPACKING=BLOCK'
    call write_block(unit,x); call write_block(unit,y); call write_block(unit,rho); call write_block(unit,ux)
    call write_block(unit,uy); call write_block(unit,temp); call write_block(unit,press)
    call write_block(unit,2.0_dp*qx); call write_block(unit,2.0_dp*qy)
    call write_block(unit,thetax); call write_block(unit,thetay); call write_block(unit,thetaz)
    call write_block(unit,sigxx); call write_block(unit,sigyy); call write_block(unit,sigxy)
    call write_block(unit,m3x); call write_block(unit,m3y); call write_block(unit,sx); call write_block(unit,sy)
    call write_block(unit,m4x); call write_block(unit,m4y); call write_block(unit,kx); call write_block(unit,ky)
    close(unit)

    call write_reconstructed_full_outputs(prefix)

    speed = sqrt(ux*ux + uy*uy)
    call choose_samples(speed, is, js, names)
    call write_vdf_samples(trim(prefix)//'_vdf_samples.dat', trim(prefix)//'_vdf_samples.txt', is, js, names)
  end subroutine write_all_outputs

  subroutine write_reconstructed_full_outputs(prefix)
    character(len=*), intent(in) :: prefix
    integer :: unit_dat, unit_csv, ii, jj, iq, jq
    real(dp) :: xf, yf, px, py, uf, vf, qxf, qyf, sigxyf
    open(newunit=unit_dat, file=trim(prefix)//'_full.dat', status='replace', action='write')
    write(unit_dat,*) 'VARIABLES = X, Y, RHO, U, V, T, P, QX_PAPER, QY_PAPER, SIGXX, SIGYY, SIGXY'
    write(unit_dat,*) 'ZONE I = ', 2*nx, ', J = ', 2*ny, ', DATAPACKING=POINT'
    open(newunit=unit_csv, file=trim(prefix)//'_full.csv', status='replace', action='write')
    write(unit_csv,'(A)') 'x,y,rho,u,v,T,p,qx_paper,qy_paper,sigxx,sigyy,sigxy'
    do jj=1,2*ny
       if (jj <= ny) then
          jq = ny-jj+1; py = -1.0_dp
       else
          jq = jj-ny; py = 1.0_dp
       end if
       yf = -0.5_dp + (real(jj,dp)-0.5_dp)/real(2*ny,dp)
       do ii=1,2*nx
          if (ii <= nx) then
             iq = nx-ii+1; px = -1.0_dp
          else
             iq = ii-nx; px = 1.0_dp
          end if
          xf = -0.5_dp + (real(ii,dp)-0.5_dp)/real(2*nx,dp)
          uf = px*ux(iq,jq)
          vf = py*uy(iq,jq)
          qxf = px*2.0_dp*qx(iq,jq)
          qyf = py*2.0_dp*qy(iq,jq)
          sigxyf = px*py*sigxy(iq,jq)
          write(unit_dat,'(12(ES23.16,1X))') xf,yf,rho(iq,jq),uf,vf,temp(iq,jq), &
               press(iq,jq),qxf,qyf,sigxx(iq,jq),sigyy(iq,jq),sigxyf
          write(unit_csv,'(*(g0,:,","))') xf,yf,rho(iq,jq),uf,vf,temp(iq,jq), &
               press(iq,jq),qxf,qyf,sigxx(iq,jq),sigyy(iq,jq),sigxyf
       end do
    end do
    close(unit_dat)
    close(unit_csv)
  end subroutine write_reconstructed_full_outputs

  subroutine choose_samples(speed, is, js, names)
    real(dp), intent(in) :: speed(nx,ny)
    integer, intent(out) :: is(6), js(6)
    character(len=24), intent(out) :: names(6)
    names = (/ 'center                  ', 'top_mid                 ', 'bottom_mid              ', 'left_mid                ', 'right_mid               ', 'max_speed               ' /)
    call nearest_cell(0.0_dp, 0.0_dp, is(1), js(1))
    call nearest_cell(0.25_dp, 0.475_dp, is(2), js(2))
    call nearest_cell(0.25_dp, 0.025_dp, is(3), js(3))
    call nearest_cell(0.025_dp, 0.25_dp, is(4), js(4))
    call nearest_cell(0.475_dp, 0.25_dp, is(5), js(5))
    call max_speed_cell(speed, is(6), js(6))
  end subroutine choose_samples

  subroutine nearest_cell(x0,y0,ii,jj)
    real(dp), intent(in) :: x0,y0
    integer, intent(out) :: ii,jj
    integer :: i,j
    real(dp) :: d, dbest
    dbest = huge(1.0_dp); ii=1; jj=1
    do j=1,ny; do i=1,nx
       d = (x(i,j)-x0)**2 + (y(i,j)-y0)**2
       if (d < dbest) then
          dbest=d; ii=i; jj=j
       end if
    end do; end do
  end subroutine nearest_cell

  subroutine max_speed_cell(speed, ii,jj)
    real(dp), intent(in) :: speed(nx,ny)
    integer, intent(out) :: ii,jj
    integer :: i,j
    real(dp) :: best
    best = -1.0_dp; ii=1; jj=1
    do j=1,ny; do i=1,nx
       if (speed(i,j) > best) then
          best = speed(i,j); ii=i; jj=j
       end if
    end do; end do
  end subroutine max_speed_cell

  subroutine write_vdf_samples(datname, txtname, is, js, names)
    character(len=*), intent(in) :: datname, txtname
    integer, intent(in) :: is(6), js(6)
    character(len=24), intent(in) :: names(6)
    integer :: unit, k, ku, kv
    real(dp) :: Ugrid(nv,nv), Vgrid(nv,nv)
    do kv=1,nv; do ku=1,nv
       Ugrid(ku,kv)=cu(ku); Vgrid(ku,kv)=cv(kv)
    end do; end do
    open(newunit=unit, file=datname, status='replace', action='write')
    write(unit,*) 'VARIABLES = U, V, H_center, H_top_mid, H_bottom_mid, H_left_mid, H_right_mid, H_max_speed, B_center, B_top_mid, B_bottom_mid, B_left_mid, B_right_mid, B_max_speed'
    write(unit,*) 'ZONE I = ', nv, ', J = ', nv, ', DATAPACKING=BLOCK'
    call write_vblock(unit,Ugrid); call write_vblock(unit,Vgrid)
    do k=1,6
       call write_vblock(unit,h(is(k),js(k),:,:))
    end do
    do k=1,6
       call write_vblock(unit,b(is(k),js(k),:,:))
    end do
    close(unit)

    open(newunit=unit, file=txtname, status='replace', action='write')
    write(unit,'(A)') 'name i j x y rho u v T'
    do k=1,6
       write(unit,'(A24,2(1X,I6),7(1X,ES23.16))') trim(names(k)), is(k), js(k), x(is(k),js(k)), y(is(k),js(k)), rho(is(k),js(k)), ux(is(k),js(k)), uy(is(k),js(k)), temp(is(k),js(k))
    end do
    close(unit)
  end subroutine write_vdf_samples

  subroutine compute_boundary_flux_diagnostic(maxflux)
    real(dp), intent(out) :: maxflux
    integer :: i,j,ku,kv
    real(dp) :: rhow, flux, fw
    maxflux = 0.0_dp
    do j=1,ny
       rhow = diffuse_rho_wall(h(nx,j,:,:),'right',0.0_dp,0.0_dp,t_cold)
       flux = 0.0_dp
       do kv=1,nv; do ku=1,nv
          if (cu(ku) > 0.0_dp) then
             fw = h(nx,j,ku,kv)
          else if (cu(ku) < 0.0_dp) then
             fw = rhow*wall_unit_M(ku,kv,0.0_dp,0.0_dp,t_cold)
          else
             fw = 0.0_dp
          end if
          flux = flux + cu(ku)*fw*w(ku,kv)
       end do; end do
       maxflux = max(maxflux,abs(flux))
    end do
    do i=1,nx
       rhow = diffuse_rho_wall(h(i,ny,:,:),'top',0.0_dp,0.0_dp,t_hot)
       flux = 0.0_dp
       do kv=1,nv; do ku=1,nv
          if (cv(kv) > 0.0_dp) then
             fw = h(i,ny,ku,kv)
          else if (cv(kv) < 0.0_dp) then
             fw = rhow*wall_unit_M(ku,kv,0.0_dp,0.0_dp,t_hot)
          else
             fw = 0.0_dp
          end if
          flux = flux + cv(kv)*fw*w(ku,kv)
       end do; end do
       maxflux = max(maxflux,abs(flux))
    end do
  end subroutine compute_boundary_flux_diagnostic

  subroutine write_metadata(iter, tnow, conv, res)
    integer, intent(in) :: iter
    real(dp), intent(in) :: tnow, res
    logical, intent(in) :: conv
    integer :: unit
    real(dp) :: maxspeed, kinetic_energy, tail_relative, rt, effective_pr
    maxspeed = maxval(sqrt(ux*ux+uy*uy))
    kinetic_energy = 4.0_dp*sum(rho*(ux*ux+uy*uy))*dx*dy
    tail_relative = exp(-min(vmin*vmin,vmax*vmax)/t_hot)
    rt = t_cold/t_hot
    if (trim(model) == 'bgk') then
       effective_pr = 1.0_dp
    else
       effective_pr = pr
    end if
    open(newunit=unit, file=trim(out_prefix)//'_metadata.json', status='replace', action='write')
    write(unit,'(A)') '{'
    write(unit,'(A,A,A)') '  "solver_version": "', solver_version, '",'
    write(unit,'(A)') '  "problem": "stationary thermal square cavity",'
    write(unit,'(A)') '  "no_external_force": true,'
    write(unit,'(A)') '  "quarter_domain": true,'
    write(unit,'(A)') '  "symmetry_planes": ["x=0 specular", "y=0 specular"],'
    write(unit,'(A,I0,A)') '  "quarter_nx": ', nx, ','
    write(unit,'(A,I0,A)') '  "quarter_ny": ', ny, ','
    write(unit,'(A,I0,A)') '  "full_reconstructed_nx": ', 2*nx, ','
    write(unit,'(A,I0,A)') '  "full_reconstructed_ny": ', 2*ny, ','
    write(unit,'(A,I0,A)') '  "nv": ', nv, ','
    write(unit,'(A,ES23.16,A)') '  "vmin": ', vmin, ','
    write(unit,'(A,ES23.16,A)') '  "vmax": ', vmax, ','
    write(unit,'(A,ES23.16,A)') '  "Kn_paper": ', kn, ','
    write(unit,'(A)') '  "kn_definition": "lambda/L = 1/(sqrt(2)*n0*pi*d^2*L)",'
    write(unit,'(A)') '  "collision_frequency_in_paper_equation": "nu=rho*sqrt(pi*T)/2",'
    write(unit,'(A)') '  "implemented_relaxation_time": "tau=2*Kn/(rho*sqrt(pi*T))",'
    write(unit,'(A,ES23.16,A)') '  "T_hot": ', t_hot, ','
    write(unit,'(A,ES23.16,A)') '  "T_cold": ', t_cold, ','
    write(unit,'(A,ES23.16,A)') '  "RT": ', rt, ','
    write(unit,'(A,ES23.16,A)') '  "Pr": ', pr, ','
    write(unit,'(A,ES23.16,A)') '  "effective_model_Pr": ', effective_pr, ','
    write(unit,'(A,ES23.16,A)') '  "dt": ', dt, ','
    write(unit,'(A,I0,A)') '  "iteration": ', iter, ','
    write(unit,'(A,ES23.16,A)') '  "time": ', tnow, ','
    write(unit,'(A,ES23.16,A)') '  "final_residual": ', res, ','
    write(unit,'(A,ES23.16,A)') '  "initial_mass_full_domain": ', initial_mass, ','
    write(unit,'(A,ES23.16,A)') '  "final_mass_full_domain": ', final_mass, ','
    write(unit,'(A,ES23.16,A)') '  "mass_drift_relative": ', mass_drift_relative, ','
    write(unit,'(A,ES23.16,A)') '  "max_physical_wall_mass_flux": ', max_physical_wall_mass_flux, ','
    write(unit,'(A,ES23.16,A)') '  "max_discrete_target_conservation_error": ', max_target_conservation_error, ','
    write(unit,'(A)') '  "target_projection": "positive exponential conservative projection",'
    write(unit,'(A,I0,A)') '  "projection_max_iter": ', projection_max_iter, ','
    write(unit,'(A,ES23.16,A)') '  "projection_tolerance": ', projection_tol, ','
    write(unit,'(A,I0,A)') '  "positive_projection_failures": ', positive_projection_failures, ','
    write(unit,'(A,I0,A)') '  "max_projection_iterations_used": ', max_projection_iterations_used, ','
    write(unit,'(A,I0,A)') '  "negative_target_values": ', negative_target_values, ','
    write(unit,'(A,I0,A)') '  "clipped_updated_values": ', clipped_updated_values, ','
    if (trim(model) == 'shakhov') then
       write(unit,'(A,ES23.16,A,ES23.16,A)') '  "shakhov_weight_limiter": [', &
            shakhov_weight_min, ',', shakhov_weight_max, '],'
       write(unit,'(A,I0,A)') '  "shakhov_limited_target_values": ', &
            shakhov_limited_target_values, ','
       write(unit,'(A,ES23.16,A)') '  "shakhov_raw_weight_min": ', shakhov_raw_weight_min, ','
       write(unit,'(A,ES23.16,A)') '  "shakhov_raw_weight_max": ', shakhov_raw_weight_max, ','
       write(unit,'(A,ES23.16,A)') '  "shakhov_limited_H_base_mass_fraction": ', &
            shakhov_limited_h_base_mass/(shakhov_total_h_base_mass+eps), ','
       write(unit,'(A,ES23.16,A)') '  "shakhov_limited_B_base_content_fraction": ', &
            shakhov_limited_b_base_content/(shakhov_total_b_base_content+eps), ','
    else
       write(unit,'(A)') '  "shakhov_weight_limiter": null,'
       write(unit,'(A)') '  "shakhov_limited_target_values": null,'
       write(unit,'(A)') '  "shakhov_raw_weight_min": null,'
       write(unit,'(A)') '  "shakhov_raw_weight_max": null,'
       write(unit,'(A)') '  "shakhov_limited_H_base_mass_fraction": null,'
       write(unit,'(A)') '  "shakhov_limited_B_base_content_fraction": null,'
    end if
    write(unit,'(A,ES23.16,A)') '  "velocity_box_maxwellian_tail_relative_at_Thot": ', tail_relative, ','
    write(unit,'(A,ES23.16,A)') '  "max_speed_nondimensional": ', maxspeed, ','
    write(unit,'(A,ES23.16,A)') '  "kinetic_energy_nondimensional": ', kinetic_energy, ','
    write(unit,'(A)') '  "heat_flux_output_convention": "paper q = integral(c^2*c*f); internal q_half is multiplied by 2 on output",'
    write(unit,'(A)') '  "quantitative_fields_are_unfiltered": true,'
    write(unit,'(A,A,A)') '  "model": "', trim(model), '",'
    if (conv) then
       write(unit,'(A)') '  "converged": true'
    else
       write(unit,'(A)') '  "converged": false'
    end if
    write(unit,'(A)') '}'
    close(unit)
  end subroutine write_metadata

  function int_to_str(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    write(s,'(I0)') i
  end function int_to_str

end program dvm2d_thermal_cavity
