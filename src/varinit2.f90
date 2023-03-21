module diheds
    contains subroutine dihedr(i1,i2,i3,i4)
    include 'mols.par'
    common/crda/x(maxatm,8)
    acdc=((180.0*7.0)/22.0)
    one=1.d0
    x1=x(i2,1)-x(i1,1)
    y1=x(i2,2)-x(i1,2)
    z1=x(i2,3)-x(i1,3)
    x2=x(i3,1)-x(i2,1)
    y2=x(i3,2)-x(i2,2)
    z2=x(i3,3)-x(i2,3)
    ux1=y1*z2-z1*y2
    uy1=z1*x2-x1*z2
    uz1=x1*y2-y1*x2
    x1=x(i4,1)-x(i3,1)
    y1=x(i4,2)-x(i3,2)
    z1=x(i4,3)-x(i3,3)
    ux2=z1*y2-y1*z2
    uy2=x1*z2-z1*x2
    uz2=y1*x2-x1*y2

    u1=ux1*ux1+uy1*uy1+uz1*uz1
    u2=ux2*ux2+uy2*uy2+uz2*uz2
    u=u1*u2

    if (u.ne.zero) then
    a=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u)
    a=max(a,-one)
    a=min(a,one)
    dihedr=acos(a)*acdc
    if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+&
       uz1*(ux2*y2-uy2*x2).lt.zero) dihedr =-dihedr
    return
    else
    write (*,'(a,4i5)')' dihedr> Error in coordinates of atoms #: '&
                      ,i1,i2,i3,i4
    stop
    endif
    end subroutine dihedr
end module diheds