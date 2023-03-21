module pdbrw
	parameter(maxstr=4096)
	parameter(maxatm=900)

	character str1*5,str2*80,str3*80,atom*4,endmdl*6
	character*3 atnam(maxstr,maxatm)
	character*4 resnam(maxstr,maxatm)
	character*1 ch(maxstr,maxatm),atn(maxstr,maxatm)
	real mene(maxstr),x(maxstr,maxatm,3),xo(maxstr,maxatm),&
        yo(maxstr,maxatm)
	integer atno(maxstr,maxatm),iresno(maxstr,maxatm),&
        imodel(maxstr)
        character iff1*200,off1*200,path*150,mname*20,NOA*17,ipdb*200,off2*200,&
        off0*200,natoms*6

contains

    subroutine readpdb

10	format(a5,4x,i4,3x,f10.2)
20	format(a4,4x,i3,2x,a2,2x,a4,i5,4x,3f8.3,a25)
30	format(a3)
40	format(a6)
50      format(a17,i5)
60	format(a4,4x,i3,2x,a3,1x,a4,i5,4x,3f8.3)
70	format(a4,4x,i3,2x,a3,1x,a4,i5,4x,3f8.3)
        
	open(unit=11,file='usermols.inp',status='old')
	read(11,'(A)')path
	read(11,*)mname
	read(11,*)
	read(11,*)ns
        close(unit=11)
        open(unit=12,file='sim.log',status='old')
        read(12,*)natoms,natom
        close(unit=12)
!    ipdb = path(:LNBLNK(path))//'/'//mname(:LNBLNK(mname))//'_'//'inf.log'
    iff1 = path(:LNBLNK(path))//'/'//mname(:LNBLNK(mname))//'_'//'tmp_mini.pdb'
    off0 = path(:LNBLNK(path))//'/'//mname(:LNBLNK(mname))//'_'//'mini.pdb'
    off1 = path(:LNBLNK(path))//'/'//mname(:LNBLNK(mname))//'_'//'minicl.pdb'
    off2 = path(:LNBLNK(path))//'/'//mname(:LNBLNK(mname))//'_'//'inpp.pdb'
	
        open(unit=1,file=iff1,status='old')

	do ij = 1,ns

	read(1,*) str1,im,ene
	read(1,*) str2

	imodel(ij) = im
	mene(im) = ene
	
	do i = 1,natom
	read(1,*) atom,atno(im,i),atnam(im,i),resnam(im,i),&
       ch(im,i),iresno(im,i),x(im,i,1),x(im,i,2),x(im,i,3),&
       xo(im,i),yo(im,i),atn(im,i)
	enddo
	
	read(1,*) str3
	read(1,*) endmdl

	enddo
	close(1)

	open(unit=2,file=off1,status='unknown')

	do i = 1,ns
	imo = imodel(i)
        open(unit=119,file='sim.log',status='old',access='append')
                write(119,*)'MODEL NO MAP:',i,imo
        close(unit=119)

!        write(2,10) str1,imo,mene(imo)
        write(2,10) str1,i,mene(imo)
        write(*,*)'MODEL ID:',(i,imo),'ENERGY (KJ):',mene(imo)
	do j = 1,natom
	write(2,20) atom,atno(imo,j),atnam(imo,j),resnam(imo,j),&
       iresno(imo,j),x(imo,j,1),x(imo,j,2),x(imo,j,3)
	enddo
	
	write(2,30)'TER'
	write(2,40)'ENDMDL'
	enddo
	close(2)


	open(unit=3,file=off2,status='unknown')

        do i = 1,1
        
        do j = 1,natom
        write(3,70) atom,atno(imo,j),atnam(imo,j),resnam(imo,j),&
        iresno(imo,j),x(imo,j,1),x(imo,j,2),x(imo,j,3)
        enddo

        enddo
        close(3)


	open(unit=6,file=off0,status='unknown')

	do i = 1,ns
	imo = imodel(i)
!	write(6,10) str1,imo,mene(imo)
	write(6,10) str1,i,mene(imo)
	
	do j = 1,natom
	write(6,60) atom,atno(imo,j),atnam(imo,j),resnam(imo,j),&
        iresno(imo,j),x(imo,j,1),x(imo,j,2),x(imo,j,3)
	enddo
	
	write(6,30)'TER'
	write(6,40)'ENDMDL'
	enddo
	close(6)
	end subroutine readpdb
end module pdbrw
