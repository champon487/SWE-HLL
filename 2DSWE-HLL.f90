
! 特性波速度の計算
subroutine propagation(sl,sc,sr,u,v,h,g,snm,spec,dd,tsc,poro,hmin)
    implicit none
    double precision, intent(in)  :: u, v, h, g, snm, spec, dd, poro, tsc, hmin
    double precision, intent(out) :: sl, sc, sr

    double precision :: cc, cf, vv, tt, aa
    
    cc = dsqrt(g*h)
    vv = dsqrt(u**2.d0+v**2.d0)

    if( vv>1e-5 .and. h>hmin ) then
        cf = g*snm**2.d0/h**(1.d0/3.d0)
        tt = cf*vv**2.d0/(spec*g*dd)
        
        if( tt>tsc ) then
            aa = 4.d0/(1.d0-poro)*(tt-tsc)**0.5d0*cc**2.d0*(1.5d0*cf*dsqrt(dd/spec/g)*u/h*(3.d0*u**2.d0+v**2.d0)/vv   &
                +(tt-tsc)*dsqrt(spec*g*dd**3.d0)*u*v**2.d0/h/(u**2.d0+v**2.d0)**1.5d0)
        else
            aa = 0.d0
        end if

    else
        cf = 0.d0
        tt = 0.d0
        aa = 0.d0
    end if

    ! Goutiere et al. (2008) ASCEによる三次方程式の近似解法

    if( u>=0.d0 ) then
        sr = u+cc
        sc = 0.5d0*(u-cc+dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
        sl = 0.5d0*(u-cc-dsqrt((u-cc)**2.d0+4.d0*aa/(u+cc)))
    else
        sr = 0.5d0*(u+cc+dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
        sc = 0.5d0*(u+cc-dsqrt((u+cc)**2.d0+4.d0*aa/(u-cc)))
        sl = u-cc
    end if

end subroutine

! HLL法によるフラックスの計算（x方向)

subroutine hllx1storder(f,ff,sl,sr,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sr, u
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j

    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                f(i,j) = (ff(i,j)*sr(i,j)-ff(i+1,j)*sl(i,j)+sr(i,j)*sl(i,j)*(u(i+1,j)-u(i,j)))/(sr(i,j)-sl(i,j))
            end if
        end do
    end do

end subroutine

! HLL法による流砂フラックスの計算（x方向)

subroutine hllzx1storder(f,ff,sl,sc,sr,z,u,snm,dx,qx,h,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: snm, dx
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, z, qx, h
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust, hst, s0dx, zzz

    do j=1,ny-1
        do i=0,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i+1,j)
            else
                ust = (u(i,j)+u(i+1,j))*0.5d0
                hst = (h(i,j)+h(i+1,j))*0.5d0
                s0dx = snm**2.*ust**2./hst**(4./3.)*dx
                zzz = (z(i+1,j)-z(i,j))+s0dx
                if( ust>0.d0 ) then
!                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i+1,j)*sl(i,j)+sc(i,j)*sl(i,j)*zzz)/(sc(i,j)-sl(i,j))
                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i+1,j)*sl(i,j)+sc(i,j)*sl(i,j)*(z(i+1,j)-z(i,j)))/(sc(i,j)-sl(i,j))
                else
!                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i+1,j)*sc(i,j)-sc(i,j)*sr(i,j)*zzz)/(sc(i,j)-sr(i,j))
                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i+1,j)*sc(i,j)-sc(i,j)*sr(i,j)*(z(i+1,j)-z(i,j)))/(sc(i,j)-sr(i,j))
                end if
            end if
        end do
    end do

end subroutine

! HLL法によるフラックスの計算(y方向)

subroutine hlly1storder(f,ff,sl,sr,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sr, u
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j

    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                f(i,j) = (ff(i,j)*sr(i,j)-ff(i,j+1)*sl(i,j)+sr(i,j)*sl(i,j)*(u(i,j+1)-u(i,j)))/(sr(i,j)-sl(i,j))
            end if
        end do
    end do

end subroutine

! HLL法による流砂フラックスの計算(y方向)

subroutine hllzy1storder(f,ff,sl,sc,sr,z,u,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, dimension(0:nx,0:ny), intent(in)  :: ff, sl, sc, sr, u, z
    double precision, dimension(0:nx,0:ny), intent(out) :: f

    integer :: i, j
    double precision :: ust

    do j=0,ny-1
        do i=1,nx-1
            if( sl(i,j)>0.d0 ) then
                f(i,j) = ff(i,j)
            else if( sr(i,j)<0.d0 ) then
                f(i,j) = ff(i,j+1)
            else
                ust = (u(i,j)+u(i,j+1))*0.5d0
                if( ust>0.d0 ) then
                    f(i,j) = (ff(i,j)*sc(i,j)-ff(i,j+1)*sl(i,j)+sc(i,j)*sl(i,j)*(z(i,j+1)-z(i,j)))/(sc(i,j)-sl(i,j))
                else
                    f(i,j) = (-ff(i,j)*sr(i,j)+ff(i,j+1)*sc(i,j)-sc(i,j)*sr(i,j)*(z(i,j+1)-z(i,j)))/(sc(i,j)-sr(i,j))
                end if
            end if
        end do
    end do

end subroutine

! 境界条件

subroutine boundary(h,qx,qy,x,y,z,dz,dis,wid,h0,snm,ib,dy,nx,ny)
    implicit none
    integer, intent(in)  :: nx, ny
    double precision, intent(in) :: dis, wid, dy, snm, ib, h0
    double precision, dimension(0:nx,0:ny), intent(in)    :: x, y
    double precision, dimension(0:nx,0:ny), intent(inout) :: h, qx, qy, z, dz

    integer :: i, j
    
    do j=0,ny
        h( 0,j) = h(   1,j)
        h(nx,j) = h(nx-1,j)

        qx( 0,j) = -qx(   1,j)     ! 上流端不透過条件
!        qx(nx,j) = -qx(nx-1,j)

!        qx( 0,j) = dis/wid

        ! 水路中央部からの流入．水深と単位幅流量の指定

        if(y(0,j)>=0.5 .and. y(0,j)<=1.d0) then
            qx(0,j) = dis/(0.5d0)
            h(0,j) = h0
        end if
        qx(nx,j) = qx(nx-1,j)

        qy( 0,j) = qy(   1,j)
        qy(nx,j) = qy(nx-1,j)

        dz( 0,j) = 0.d0
        dz(nx,j) = dz(nx-1,j)

        z( 0,j) = z( 0,j)+dz( 0,j)
        z(nx,j) = z(nx,j)+dz(nx,j)
    end do

    do i=0,nx
        h(i, 0) = h(i,   1)
        h(i,ny) = h(i,ny-1)

        qx(i, 0) = qx(i,   1)
        qx(i,ny) = qx(i,ny-1)

        qy(i, 0) = -qy(i,   1)
        qy(i,ny) = -qy(i,ny-1)

        dz(i, 0) = dz(i,   1)
        dz(i,ny) = dz(i,ny-1)

        z(i, 0) = z(i,   1)
        z(i,ny) = z(i,ny-1)
    end do

end subroutine

! Paraview可視化用のファイル出力

subroutine out2paraview(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: tt
    double precision, dimension(0:nx,0:ny), intent(in)  :: x, y, z, h, u, v, z0
    character(20),intent(in) :: fpout

    integer ::  i, j

    open(100,file=fpout,status='unknown')
				
		write(100,'(a26)') '# vtk DataFile Version 3.0'
		write(100,'(a13,f10.3,a7)') 'Time= ', tt
		write(100,'(a5)') 'ASCII'
		write(100,'(a23)') 'DATASET STRUCTURED_GRID'
		write(100,'(a10,3i8)') 'DIMENSIONS', nx+1, ny+1, 1
		write(100,'(a6,i15,a7)') 'POINTS', (nx+1)*(ny+1), 'double'
			
		do j=0,ny
			do i=0,nx
				write(100,*) x(i,j), y(i,j), 0.d0
			end do
		end do
				
		write(100,'(a10,i15)') 'POINT_DATA', (nx+1)*(ny+1)
		write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'H', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)+h(i,j)
			end do
		end do
        
		write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'hs', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) h(i,j)
			end do
		end do

        write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'Z', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)
			end do
		end do
        
        write(100,'(a7,2x,a5,2x,a7,i4)') 'SCALARS', 'DZ', 'double', 1
		write(100,'(a20)') 'LOOKUP_TABLE default'
			
		do j=0,ny
			do i=0,nx
				write(100,*) z(i,j)-z0(i,j)
			end do
		end do

        write(100,'(a23)') 'VECTORS velocity double'
        do j=0,ny
            do i=0,nx
                write(100,*) u(i,j),v(i,j),0.
            end do
        end do
				
	close(100) 

end subroutine

! テキスト出力

subroutine out2txt(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
    implicit none
    integer, intent(in) :: nx, ny
    double precision, intent(in) :: tt
    double precision, dimension(0:nx,0:ny), intent(in)  :: x, y, z, h, u, v, z0
    character(20),intent(in) :: fpout

    integer ::  i, j

    open(101,file=fpout,status='unknown')
				
		write(101,*) tt, nx, ny

        do j=0,ny
            do i=0,nx
                write(101,*) x(i,j), y(i,j), h(i,j)
            end do
        end do
			
	close(101) 

end subroutine

program SWE_HLL

    implicit none 
    integer :: i, j, nx, ny, m, ii
    double precision :: dx, dy, dt, rdx, rdy, ct, dt1, dt2, hmin, h0, v0
    double precision :: g, nu, dis, chlen, wid, snm, ib, spec, diam, poro, tuk, etime, bedtime, z00
    double precision :: hmax, tt, optime, ss, cl, cr, ust, cst
    double precision :: pressure, roughness, advection, volume
    double precision :: tsc, tau, vv, dzdx, dzdy, dzdn, qbs, qbn, uc, vc, hc, frc, mu_s, sigl, sigr, pi, wdz
    double precision :: wsl1, wsc1, wsr1, wsl2, wsc2, wsr2
    double precision, dimension(:,:), allocatable :: x, y, z, wz, z0, dz, h, wh, u, v, qx, wqx, qy, wqy, sr, sl, sc
    double precision, dimension(:,:), allocatable :: f1, f2, f3, ffx, ffy, qbsta
    double precision, dimension(:,:), allocatable :: qbx, qby, qbx1, qby1, qbx2, qby2, f4
    double precision, dimension(:,:), allocatable :: qsta, fsta, dh
    character(20) :: fpout, fptxt

    g  = 9.81d0
    nu = 1e-6
    pi = 3.14159d0

    h0 = 0.096d0                        ! 上流端水深 m
    v0 = 1.94d0                         ! 上流端流速 m/s

    dis     = h0*v0*0.5d0               ! 流量 m3/s
    chlen   = 10.d0                     ! 水路長さ m
    wid     = 1.5d0                     ! 水路幅 m
    snm     = 0.01d0                    ! マニング粗度係数
    ib      = 0.0d0                     ! 水路勾配
    spec    = 1.65d0                    ! 河床材料水中比重
    diam    = 0.76d0                    ! 河床材料粒径 m
    poro    = 0.4d0                     ! 河床空隙率
    tsc     = 0.d0                      ! 限界掃流力
    mu_s    = 0.7d0                     ! 土粒子摩擦係数
    hmin    = 0.0001d0                  ! 計算上の最小水深 m

    nx = 800                            ! 計算格子数 nx→x方向，ny→y方向
    ny = 120
    dx = chlen/dble(nx)
    dy = wid/dble(ny)
    dt = 0.0001d0
    ct = 0.1d0                          ! クーラン数
    
    rdx = 1.d0/dx
    rdy = 1.d0/dy
    
    tuk   = 0.2                         ! 計算結果出力間隔　s
    bedtime = 300.
    etime = tuk*150.

    allocate( x(0:nx,0:ny), y(0:nx,0:ny), z(0:nx,0:ny), wz(0:nx,0:ny), z0(0:nx,0:ny), dz(0:nx,0:ny), h(0:nx,0:ny), wh(0:nx,0:ny), qbsta(0:nx,0:ny) )
    allocate( u(0:nx,0:ny), v(0:nx,0:ny), qx(0:nx,0:ny), wqx(0:nx,0:ny), qy(0:nx,0:ny), wqy(0:nx,0:ny) )
    allocate( sr(0:nx,0:ny), sl(0:nx,0:ny), sc(0:nx,0:ny), ffx(0:nx,0:ny), ffy(0:nx,0:ny), f1(0:nx,0:ny), f2(0:nx,0:ny), f3(0:nx,0:ny) )
    allocate( qbx(0:nx,0:ny), qby(0:nx,0:ny), qbx1(0:nx,0:ny), qby1(0:nx,0:ny), qbx2(0:nx,0:ny), qby2(0:nx,0:ny), f4(0:nx,0:ny) )
    allocate( qsta(0:nx,0:ny), fsta(0:nx,0:ny), dh(0:nx,0:ny) )

    ! 座標と河床高さの設定

    do j=0,ny
        x(0,j) = 0.d0
        do i=1,nx
            x(i,j) = x(i-1,j)+dx
        end do
    end do
    
    do i=0,nx
        y(i,0) = 0.d0
        do j=1,ny
            y(i,j) = y(i,j-1)+dy
        end do
    end do
       
    z00  = chlen*ib

    do j=0,ny
        z(0,j) = z00
        do i=1,nx
            z(i,j) = z(i-1,j)-ib*dx
            h(i,j) = hmin
            qx(i,j) = 0.d0
            qy(i,j) = 0.d0
        end do
    end do

    dz = 0.d0

    call boundary(h,qx,qy,x,y,z,dz,dis,wid,h0,snm,ib,dy,nx,ny)

    do j=0,ny
        do i=0,nx
            if( h(i,j)>hmin ) then
                u(i,j) = qx(i,j)/h(i,j)
                v(i,j) = qy(i,j)/h(i,j)
            else
                u(i,j) = 0.d0
                v(i,j) = 0.d0
            end if
        end do
    end do

    wh  = h
    wz  = z
    wqx = qx
    wqy = qy 
    
    tt = 0.
    optime = 0.
    
	m = 0
    fpout = 'vc000.vtk'
    fptxt = 'tt000.txt'

    ! 時間方向ループ

    do  
        ! 第一段階　d/dxの項
        
        h  = wh
        z  = wz
        qx = wqx
        qy = wqy

        ! 格子点のフラックスなど計算
        
        do j=0,ny
            do i=0,nx
                if ( wh(i,j)>hmin ) then
                    u(i,j) = wqx(i,j)/wh(i,j)
                    v(i,j) = wqy(i,j)/wh(i,j)
                    
                    ffx(i,j) = wqx(i,j)**2.d0/wh(i,j)+0.5d0*g*wh(i,j)**2.d0     !+g*wh(i,j)*z(i,j)
                    ffy(i,j) = wqx(i,j)*wqy(i,j)/wh(i,j)

                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*wh(i,j)**(1.d0/3.d0))

                    if( tau>tsc ) then
                        qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                        qbx1(i,j) = u(i,j)/vv*qbs
                    else
                        qbx1(i,j) = 0.d0
                    end if
                else
                    u(i,j) = 0.d0
                    v(i,j) = 0.d0
                    ffx(i,j) = 0.d0
                    ffy(i,j) = 0.d0
                    qbx1(i,j) = 0.d0
                end if

            end do
        end do

        ! 格子境界の計算
        
        do j=1,ny-1
            do i=0,nx-1
                uc = (u(i,j)+u(i+1,j))*0.5d0
                vc = (v(i,j)+v(i+1,j))*0.5d0
                hc = (wh(i,j)+wh(i+1,j))*0.5d0

                vv = dsqrt(uc**2.d0+vc**2.d0)

                if ( hc>hmin .and. vv>1e-5 ) then
                    tau = snm**2.d0*vv**2.d0/(spec*diam*hc**(1.d0/3.d0))
                else
                    tau = 0.d0
                end if

                dzdx = (-z(i,j)+z(i+1,j))*rdx
                dzdy = (-(z(i,j-1)+z(i+1,j-1))+(z(i,j+1)+z(i+1,j+1)))*0.25d0*rdy
                dzdn = (-vc*dzdx+uc*dzdy)/vv

                if( tau>tsc ) then
                    qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
                    qbx2(i,j) = -vc/vv*qbn
                else
                    qbx2(i,j) = 0.d0
                end if

            end do
        end do

        ! 特性速度の計算

        do j=1,ny-1
            do i=0,nx-1
                call propagation(wsl1, wsc1, wsr1,u(i,j),v(i,j),wh(i,j),g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl2, wsc2, wsr2,u(i+1,j),v(i+1,j),wh(i+1,j),g,snm,spec,diam,tsc,poro,hmin)
                sl(i,j) = wsl1  !(wsl1+wsl2)*0.5d0
                sc(i,j) = (wsc1+wsc2)*0.5d0
                sr(i,j) = wsr2  !(wsr1+wsr2)*0.5d0
            end do
        end do

        dt1 = 999.d0
        do j=1,ny-1
            do i=0,nx-1
                dt1 = min(dabs(dx/sl(i,j)), dabs(dx/sc(i,j)), dabs(dx/sr(i,j)),dt1)
            end do
        end do

        ! 連続式，x方向運動方程式のフラックス計算

        call hllx1storder(f1,wqx,sl,sr,wh,nx,ny)
        call hllx1storder(f2,ffx,sl,sr,wqx,nx,ny)
        
        ! y方向運動方程式のフラックス計算（ここはdoner cellのような計算．HLLは使わない)

        do j=1,ny-1
            do i=0,nx-1
                hc = (wh(i,j)+wh(i+1,j))*0.5d0
                if( hc>hmin ) then
                    ust = f1(i,j)/hc
                    if( ust>0 ) then
                        f3(i,j) = ust*wqy(i,j)
                    else
                        f3(i,j) = ust*wqy(i+1,j)
                    end if
                else
                    f3(i,j) = 0.d0
                end if

            end do
        end do

        ! 流砂フラックスの計算（固定床計算では使わない）
        
        call hllzx1storder(f4,qbx1,sl,sc,sr,wz,u,snm,dx,wqx,wh,nx,ny)
    
        ! 連続式計算

        do j=1,ny-1
            do i=1,nx-1
                h(i,j) = wh(i,j)-(-f1(i-1,j)+f1(i,j))*rdx*dt
                if ( h(i,j)<=hmin ) h(i,j) = hmin
            end do
        end do

        ! 運動方程式計算
            
        do j=1,ny-1
            do i=1,nx-1
                pressure = -g*wh(i,j)*(-z(i-1,j)+z(i+1,j))*rdx*0.5d0
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/wh(i,j)**(4.d0/3.d0)

                sigl = f2(i-1,j)-(sr(i-1,j))/(sr(i-1,j)-sl(i-1,j))*g*(wh(i,j)+wh(i-1,j))*0.5*(wz(i,j)-wz(i-1,j))
                sigr = f2(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(wh(i+1,j)+wh(i,j))*0.5*(wz(i+1,j)-wz(i,j))

                advection = (-sigl+sigr)*rdx

                qx(i,j) = (wqx(i,j)+(-advection)*dt)/(1.0+roughness*dt)
            end do
        end do

        do j=1,ny-1
            do i=1,nx-1
                qy(i,j) = wqy(i,j)-(-f3(i-1,j)+f3(i,j))*rdx*dt
            end do
        end do
        
        if( tt>bedtime ) then

            do j=1,ny-1
                do i=1,nx-1
                    wdz = -((-f4(i-1,j)+f4(i,j))*rdx+(-qbx2(i-1,j)+qbx2(i,j))*rdx)*dt/(1.d0-poro)
                    z(i,j) = wz(i,j)    !+wdz       ! 固定床計算では河床変動はさせない
                    dz(i,j) = wdz
                end do
            end do

        end if

        
        ! 第二段階　d/dyの項

        call boundary(h,qx,qy,x,y,z,dz,dis,wid,h0,snm,ib,dy,nx,ny)

        ! 格子点のフラックスなど計算
        
        do j=0,ny
            do i=0,nx
                if ( h(i,j)>hmin ) then
                    u(i,j) = qx(i,j)/h(i,j)
                    v(i,j) = qy(i,j)/h(i,j)
                    
                    ffx(i,j) = qx(i,j)*qy(i,j)/h(i,j)
                    ffy(i,j) = qy(i,j)**2.d0/h(i,j)+0.5d0*g*h(i,j)**2.d0    !+g*h(i,j)*z(i,j)
    
                    vv = dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)
                    tau = snm**2.d0*vv**2.d0/(spec*diam*h(i,j)**(1.d0/3.d0))

                    if( tau>tsc ) then
                        qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                        qby1(i,j) = v(i,j)/vv*qbs
                    else
                        qby1(i,j) = 0.d0
                    end if
                else
                    u(i,j) = 0.d0
                    v(i,j) = 0.d0
                    ffx(i,j) = 0.d0
                    ffy(i,j) = 0.d0
                    qby1(i,j) = 0.d0
                end if
            end do
        end do
        
        ! 格子境界の計算
                
        do j=0,ny-1
            do i=1,nx-1
                uc = (u(i,j)+u(i,j+1))*0.5d0
                vc = (v(i,j)+v(i,j+1))*0.5d0
                hc = (h(i,j)+h(i,j+1))*0.5d0
                
                vv = dsqrt(uc**2.d0+vc**2.d0)

                if ( hc>hmin .and. vv>1e-5 ) then
                    tau = snm**2.d0*vv**2.d0/(spec*diam*hc**(1.d0/3.d0))
                else
                    tau = 0.d0
                end if

                dzdx = (-(z(i-1,j)+z(i-1,j+1))+(z(i+1,j)+z(i+1,j+1)))*0.25d0*rdx
                dzdy = (-z(i,j)+z(i,j+1))*rdy
                dzdn = (-vc*dzdx+uc*dzdy)/vv

                if( tau>tsc ) then
                    qbs = 4.d0*(tau-tsc)**1.5d0*dsqrt(spec*g*diam**3.d0)
                    qbn = -qbs*dsqrt(tsc/tau)/mu_s*dzdn
                    qby2(i,j) = uc/vv*qbn
                else
                    qby2(i,j) = 0.d0
                end if
            end do
        end do

        do i=0,nx
!            qby1(i, 0) = -qby1(i,   1)
!            qby1(i,ny) = -qby1(i,ny-1)
            qby1(i, 0) = 0.d0
            qby1(i,ny) = 0.d0
        end do

        ! 特性波速度の計算

        do j=0,ny-1
            do i=1,nx-1
                call propagation(wsl1, wsc1, wsr1,v(i,j),u(i,j),h(i,j),g,snm,spec,diam,tsc,poro,hmin)
                call propagation(wsl2, wsc2, wsr2,v(i,j+1),u(i,j+1),h(i,j+1),g,snm,spec,diam,tsc,poro,hmin)
                sl(i,j) = wsl1  !(wsl1+wsl2)*0.5d0
                sc(i,j) = (wsc1+wsc2)*0.5d0
                sr(i,j) = wsr2  !(wsr1+wsr2)*0.5d0
            end do
        end do

        dt2 = 999.d0
        do j=0,ny-1
            do i=1,nx-1
                dt2 = min(dabs(dy/sl(i,j)), dabs(dy/sc(i,j)), dabs(dy/sr(i,j)),dt2)
            end do
        end do

        ! HLLによる連続式，y方向運動方程式のフラックス計算

        call hlly1storder(f1,qy,sl,sr,h,nx,ny)
        call hlly1storder(f3,ffy,sl,sr,qy,nx,ny)

        do j=0,ny-1
            do i=1,nx-1
                hc = (h(i,j)+h(i,j+1))*0.5d0
                if( hc>hmin ) then
                    ust = f1(i,j)/hc
                    if( ust>0 ) then
                        f2(i,j) = ust*qx(i,j)
                    else
                        f2(i,j) = ust*qx(i,j+1)
                    end if
                else
                    f2(i,j) = 0.d0
                end if

            end do
        end do
        
        ! 流砂フラックスの計算（固定床計算では使わない）

        call hllzy1storder(f4,qby1,sl,sc,sr,z,v,nx,ny)
        
        ! 連続式の計算

        do j=1,ny-1
            do i=1,nx-1
                wh(i,j) = h(i,j)-(-f1(i,j-1)+f1(i,j))*rdy*dt
                if ( wh(i,j)<=hmin ) wh(i,j) = hmin
            end do
        end do

        ! 運動方程式の計算
            
        do j=1,ny-1
            do i=1,nx-1
                pressure = -g*h(i,j)*(-z(i,j-1)+z(i,j+1))*rdy*0.5d0
                roughness = g*snm**2.d0*dsqrt(u(i,j)**2.d0+v(i,j)**2.d0)/h(i,j)**(4.d0/3.d0)

                sigl = f3(i,j-1)-(sr(i,j-1))/(sr(i,j-1)-sl(i,j-1))*g*(h(i,j)+h(i,j-1))*0.5*(z(i,j)-z(i,j-1))
                sigr = f3(i,j)-(sl(i,j))/(sr(i,j)-sl(i,j))*g*(h(i,j+1)+h(i,j))*0.5*(z(i,j+1)-z(i,j))

                advection = (-sigl+sigr)*rdy

                wqy(i,j) = (qy(i,j)+(-advection)*dt)/(1.0+roughness*dt)
            end do
        end do

        do j=1,ny-1
            do i=1,nx-1
                wqx(i,j) = qx(i,j)-(-f2(i,j-1)+f2(i,j))*rdy*dt
            end do
        end do
                
        if( tt>bedtime ) then

            do j=1,ny-1
                do i=1,nx-1
                    wdz = -((-f4(i,j-1)+f4(i,j))*rdy+(-qby2(i,j-1)+qby2(i,j))*rdy)*dt/(1.d0-poro)
                    wz(i,j) = z(i,j)    !+wdz       ! 固定床計算では河床変動させない
                    dz(i,j) = wdz
                end do
            end do

        end if

        call boundary(wh,wqx,wqy,x,y,wz,dz,dis,wid,h0,snm,ib,dy,nx,ny)

        ! 計算時間刻みの更新

        dt = min(dt1,dt2)*ct

        ! 計算結果の出力
            
        if (optime>tuk .or. m==0) then

            volume = 0.

            do j=1,ny-1
                do i=1,nx-1
                    volume = volume+h(i,j)
                end do
            end do

            write(*,'(f12.5,f10.6,e17.5)') tt, dt, volume

            write(fpout(3:5),'(i3)') m
            
            do ii=3,4
                if(fpout(ii:ii) == ' ') fpout(ii:ii) = '0'
            end do

            write(fptxt(3:5),'(i3)') m
            
            do ii=3,4
                if(fptxt(ii:ii) == ' ') fptxt(ii:ii) = '0'
            end do

            call out2paraview(x,y,z,h,u,v,z0,nx,ny,fpout,tt)
            call out2txt(x,y,z,h,u,v,z0,nx,ny,fptxt,tt)

            if( m>0 ) optime = optime-tuk
            m = m+1
        
        end if
        
        optime = optime+dt
        tt = tt+dt

        if( tt>etime ) exit
        
    end do

end program SWE_HLL