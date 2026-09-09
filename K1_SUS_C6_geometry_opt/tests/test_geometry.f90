program test_geometry
  use mod_types
  use mod_config
  use mod_geometry
  use mod_material
  use mod_operator, only: fast_pair=>pair_map_local, fast_visible=>visible_pair
  use mod_operator_old, only: old_pair=>pair_map_local, old_visible=>visible_pair
  use mod_operator_exact, only: exact_pair=>pair_map_local
  use mod_c6_modal, only: fast_step=>apply_c6_modal_operator_matrix_free
  use mod_modal_exact, only: exact_step=>apply_c6_modal_operator_matrix_free
  use omp_lib
  implicit none
  type(sim_config_type) :: c
  type(panel_mesh_type) :: mesh
  complex(dp) :: a(3,3),b(3,3),old(3,3)
  complex(dp), allocatable :: j(:,:,:),jf(:,:,:),jr(:,:,:)
  real(dp) :: x,y,z,t,rel,maxerr,t0,tf,tr
  integer :: variant,pt,ps,qt,qs,k,missed,changed,blocked,checked,i
  logical :: vf,vo,dense
  call load_config('config_smoke.nml',c)
  call load_material('config_smoke.nml',c)
  do variant=1,3
    c%n_face=3; c%n_z_cone=4; c%n_z_pipe=3
    if(variant==2) c%d_in=c%d_out
    if(variant==3) c%d_in=0.040_dp
    call build_c6_mesh(c,mesh)
    missed=0; changed=0; blocked=0; checked=0; maxerr=0
    do pt=1,6
      do ps=1,6
        do qt=1,mesh%Q
          do qs=1,mesh%Q
            vf=fast_visible(mesh%x(pt,qt),mesh%y(pt,qt),mesh%z(pt,qt), &
                            mesh%x(ps,qs),mesh%y(ps,qs),mesh%z(ps,qs),c)
            vo=old_visible(mesh%x(pt,qt),mesh%y(pt,qt),mesh%z(pt,qt), &
                           mesh%x(ps,qs),mesh%y(ps,qs),mesh%z(ps,qs),c)
            if(vo .and. .not.vf) missed=missed+1
            if(.not.vo .and. vf) error stop 'Exact visibility accepted old blocked segment'
            if(.not.vf) blocked=blocked+1
            ! Independently sample 101 points, including both endpoints.
            dense=.true.
            do k=0,100
              t=real(k,dp)/100
              x=(1-t)*mesh%x(ps,qs)+t*mesh%x(pt,qt)
              y=(1-t)*mesh%y(ps,qs)+t*mesh%y(pt,qt)
              z=(1-t)*mesh%z(ps,qs)+t*mesh%z(pt,qt)
              if(.not.inside_hexagon(x,y,z,c,0.0_dp)) dense=.false.
            end do
            if(vf .and. .not.dense) error stop 'Visible segment leaves cavity'
            call fast_pair(mesh,c,pt,qt,ps,qs,a)
            call exact_pair(mesh,c,pt,qt,ps,qs,b)
            if(maxval(abs(a-b))>1e-12_dp*max(1.0_dp,maxval(abs(b)))) error stop 'Culling changed a nonzero map'
            maxerr=max(maxerr,maxval(abs(a-b)))
            call old_pair(mesh,c,pt,qt,ps,qs,old)
            if(maxval(abs(a-old))>1e-12_dp*max(1.0_dp,maxval(abs(old)))) changed=changed+1
            checked=checked+1
          end do
        end do
      end do
    end do
    allocate(j(3,mesh%Q,6),jf(3,mesh%Q,6),jr(3,mesh%Q,6))
    do k=1,6
      do i=1,mesh%Q
        j(:,i,k)=[cmplx(sin(real(i*k,dp)),cos(real(i+k,dp)),dp), &
                    cmplx(cos(real(i*k,dp)),sin(real(i+k,dp)),dp),cmplx(0.3_dp,-0.7_dp,dp)]
      end do
    end do
    t0=omp_get_wtime()
    call exact_step(mesh,c,[0,1,2,3,4,5],j,jr)
    tr=omp_get_wtime()-t0
    t0=omp_get_wtime()
    call fast_step(mesh,c,[0,1,2,3,4,5],j,jf)
    tf=omp_get_wtime()-t0
    rel=sqrt(sum(abs(jf-jr)**2)/max(sum(abs(jr)**2),tiny(1.0_dp)))
    if(rel>1e-12_dp) error stop 'C6 block exclusion changed current'
    print *, 'variant,pairs,blocked,old_visibility_misses,changed_maps=',variant,checked,blocked,missed,changed
    print *, 'max_map_difference,relative_C6_difference=',maxerr,rel
    print *, 'exact_unculled_seconds,optimized_seconds=',tr,tf
    deallocate(j,jf,jr)
  end do
  ! Narrow obstruction at the taper/pipe seam: old ten-point sampling misses.
  c%d_in=0.180_dp; c%d_out=0.065_dp
  x=radius_at_z(0.0_dp,c)*cos(PI/6.0_dp)-0.005_dp
  y=radius_at_z(c%l_cone,c)*cos(PI/6.0_dp)-0.0001_dp
  vf=fast_visible(y,0.0_dp,c%l_cone+0.001_dp,x,0.0_dp,0.0_dp,c)
  vo=old_visible(y,0.0_dp,c%l_cone+0.001_dp,x,0.0_dp,0.0_dp,c)
  if(vf .or. .not.vo) error stop 'Narrow junction obstruction test failed'
  print *, 'Geometry / culling / C6 equivalence tests PASS'
end program
