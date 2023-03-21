    
    
    
    
    subroutine energy

	include 'mols.par'
        common/crda/x(maxatm,8) 
        common/ranges/jstart(maxatm,10),jend(maxatm,10),j1_4(maxatm,12)
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
        common/crdb/y(maxatm,8) 
        common/vectors/iv(maxpar,4)
        common /par/ natom, ntor, nhb, ns, lres
        common/hb/ihb1(maxhb),ihb2(maxhb),c(maxhb),d(maxhb)
        common/order/n
        common/inp/e_inp(maxstr),p1(maxstr,maxpar)
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)

	common /ctl/iopt,iff,icint,fseq(maxres)
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4

	real fret,ftol,p(maxpar)
	integer iter, n



    call getx
        do i=1,natom
          do j=4,8
             y(i,j)=x(i,j)
          enddo          
	enddo
ccccc	call pdb2
	do i = 1,ns
	  do j=1,n
	    p(j) = p1(i,j)
        enddo
    call frprmn(p,n,ftol,iter,fret)
    call outp(p,fret,i)
c	write(6,*) 'Struc. NO:	Energy value	No. of iteration'	
  write(6,*) 'minimized structure no: ',i,fret,iter
c*******Flushing standard output************
  call flush(6)
c*******End of Flushing ********************
  write(31,*) 'minimized structure no: ',i,fret,iter
  enddo
c	call pdbout
  return
  end


  subroutine getx
	include 'mols.par'
        common/inp/e_inp(maxstr),p1(maxstr,maxpar)
        common /par/ natom, ntor, nhb, ns, lres
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
        common/order/n
10	format(1x,8(1x,f5.1),5x,f20.4)
	open(unit=2,file='mols_q.out',status='unknown')
	do i = 1,ns
cc	e_inp(i) = emo(i)
	do j = 1,n
cc	p1(i,j) = opmo(i,j)
	enddo
c	print *,(p1(i,j1),j1=1,n),e_inp(i) 
cc	read(2,10)(p1(i,j),j=1,n),e_inp(i) 
	read(2,*)(p1(i,j),j=1,n),e_inp(i)
	enddo
	close(unit=2)
	return
	end
c********************************************************************
c*********************************************************************
	function func(p2)
	include 'mols.par'
        common/order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,iff,icint,fseq(maxres)
	dimension p2(maxpar)
c
        do j = 1,n
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
        enddo
c
	call molgen(p2)
c	call energee2(ej)
c	func = ej
c	func = ampene(1,1)
		if(iff.eq.1) func=ampene(1,1)
		if(iff.eq.2) func=ecpene(1,1)

	return
	end
c*******************************************************************
	subroutine dfunc(p2,g1)
	include 'mols.par'
        common/order/n
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ctl/iopt,iff,icint,fseq(maxres)	
	dimension p2(maxpar),g1(maxpar)
	real deltaX,deltaE,f1,p2
	deltaX = 0.1
c
        do j = 1,n
        p2(j) = mod(p2(j),360.0)
        if(p2(j).le.0) p2(j) = 360.0 + p2(j)
	phi(j) = p2(j)
        enddo
c
	call molgen(p2)
c	call energee2(ej)
c	f1 = ej
		if(iff.eq.1) f1=ampene(1,1)
		if(iff.eq.2) f1=ecpene(1,1)
ccccc	f1 = ecpene(1,1)
	do j=1,n
	g1(j) = 0.0
	deltaE = 0.0
	ej = 0.0
	p2(j) = p2(j) - deltaX
	phi(j) = p2(j)
	call molgen(p2)
c	call energee2(ej)
c	ej = ampene(1,1)
		if(iff.eq.1) ej=ampene(1,1)
		if(iff.eq.2) ej=ecpene(1,1)
	deltaE = f1-ej
	g1(j) = deltaE/deltaX
	p2(j) = p2(j) + deltaX
	phi(j) = p2(j)
	enddo

c214	format (/,4(2x,2hG(,i1,4h) = ,f10.4))
c	write (6,214) (j,g1(j), j=1,n)

	return
	end



    subroutine outp(p3,ene,i)
        include 'mols.par'
            common/order/n
            common /par/ natom, ntor, nhb, ns, lres
        common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
         &              emo(maxstr),emi(maxstr)
        common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
        character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
        dimension p3(maxpar)
    
    10	format(i4,1x,f8.3,<n>(1x,f6.1))
            do k1=1,n
        opmi(i,k1) = p3(k1)
        p3(k1) = p3(k1)-180.0
        enddo 
    c	print *,i,(opmi(i,k1),k1=1,n)
        emi(i) =ene
        open(unit=14,file=of1,status='unknown')
    c           write(14,*)(p3(j),j=1,n),ene 
    c           write(14,10)i,ene,(p3(j),j=1,n)
               write(14,10)i,ene,(p3(j),j=1,n)
    c10	format(i4,1x,f15.2,<n>(1x,f6.1))
        if(i.eq.ns) close(unit=14)
            return
            end
