c     library name : ecppar.f   

c     Copyright (c) 2013      K.Vengadesan and N.Gautham

c     This library is free software; you can redistribute it and/or
c     modify it under the terms of the GNU Lesser General Public
c     License as published by the Free Software Foundation; either
c     version 2.1 of the License, or (at your option) any later version.
c
c     This library is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c     Lesser General Public License for more details.
c
c     You should have received a copy of the GNU Lesser General Public
c     License along with this library; if not, write to the Free
c     Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
c     02110-1301  USA
c
c
c     contact : n_gautham@hotmail.com  



c	program to generate MOLS parameter file from the pdb file
c	program pargen
	subroutine ecppar(iopt,fseq,itcount)
	parameter (maxbnd = 815, jx = 25)
	include 'mols.par'
	include 'ecepp.h'
	common /par/ natom, ntor, nhb, ns, lres
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
c
      dimension hbc(mxhbdo,mxhbac),hba(mxhbdo,mxhbac)
      logical do(mxtyat),ac(mxtyat)
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	character atname(maxatm)*4, rename(maxatm)*4, tcode*4, scode*1,
     &  an(maxatm,jx)*4, bn(maxatm,jx)*4, a(maxatm,jx)*4, 
     &  b(maxatm,jx)*4, c(maxatm,jx)*4, cn(maxatm,jx)*4,hbat(maxatm)*1, 
     &  dn(maxatm,jx)*4, m14(maxatm,jx)*4,atype*2, atnam*2, at1*4, 
     &  at2*1, at3*2, d1*4, d2*4, d3*4, d4*4,hbat1(maxatm)*2,
     &  fseq(maxatm)*1
c
	integer iresno(maxatm), atmno(maxatm), na(maxatm,jx), ka(maxatm),
     &  kc(maxatm),ma(maxatm,jx), kb(maxatm), mb(maxatm,jx), 
     &  nb(maxatm,jx), mc(maxatm,jx), nc(maxatm,jx), kd(maxatm), 
     &  nd(maxatm,maxatm), ke(maxatm), ne(maxatm,jx),k14(maxatm),
     &  k15(maxatm), n15(maxatm,jx), md(maxatm,jx), n14(maxatm,jx),
     &  kt(maxatm),kf15(maxatm),nt(maxatm,maxatm),nf15(maxatm,jx)

	real ele
cc	write(31,*)'Enter No. of atoms'
cc	read *,natom
c
10	format(8x,i3,1x,a4,1x,a4,3x,i3)
20	format(a4,a1,1x,i2)
30	format(8x,2a4)
40	format(10i4)
50	format(i4,1x,f3.1,1x,f3.1)
60	format(2x,20i4)
70	format(2x,a2,6x,f6.3,5x,f6.3)
80	format(a4,7x,a1,2x,a2,f8.4,f7.3,2x,i2)
90	format(a4,a1,1x,i2)
100	format(i4,1x,f8.4,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1)
110	format(8x,4a4)
120	format(i4,i4,i6,i6)
c
c	read atom number, atom name, residue name and residue number
c	from the pdb file generated by pdbgen
c
	open(unit=1,file=pf1,status='old')
	do i=1,natom
		read(1,10) atmno(i),atname(i),rename(i),iresno(i)
c		write(31,*)atmno(i),atname(i),rename(i),iresno(i)
	enddo
	close(unit=1)
c****************************to find first neighbours*************************
      do i=1,natom
	ka(i)=1
	open(unit=2,file='ALLCONN.lib',status='old')
	rewind 2
	do i1=1,maxbnd
	  read(2,20)tcode,scode,iatom
	  if(tcode.eq.rename(i))then
	   do l=1,iatom
	    read(2,30) (a(i,j),j=1,2)
	    if(a(i,1).eq.atname(i))then
	      do j=1,natom
		if(atname(j).eq.a(i,2).and.a(i,2).ne.'    ')then
	 	  if(a(i,1).eq.' C  '.and.a(i,2).eq.' N  ')then
		     itemp=iresno(i)+1
		  else
	 	     itemp=iresno(i)
		  endif
	          if(iresno(j).eq.itemp)then
			ma(i,1)=i
			ma(i,2)=j
c			write(31,*)i,j
		  endif
	        endif
	      enddo
c	write(31,*) (ma(i,j1),j1=1,2), (a(i,j1),j1=1,2)
	      if(ma(i,2).ne.na(i,1).and.ma(i,2).ne.na(i,2).and.
     &	      ma(i,2).ne.na(i,3).and.ma(i,2).ne.na(i,4))then
		na(i,ka(i))=ma(i,2)
		an(i,ka(i))=a(i,2)
		ka(i)=ka(i)+1
	      endif
	    endif
	  enddo
	 endif
	enddo
c	write(6,*) i,(an(i,j1),j1=1,ka(i)-1)
c	write(6,*) i,(na(i,j1),j1=1,ka(i)-1)
c*************************to find second neighbours**************************
	kb(i)=1
	do l2=1,ka(i)-1
	rewind 2
	  do i2=1,maxbnd
	    read(2,20)tcode,scode,iatom
	    if(tcode.eq.rename(i))then
	      do l1=1,iatom
		read(2,30) (b(i,j),j=1,2)
		if(an(i,l2).eq.b(i,1))then
		  do j=1,natom
	 	    if(atname(j).eq.b(i,2).and.b(i,2).ne.'    ')then
		      if(b(i,1).eq.' C  '.and.b(i,2).eq.' N  ')then
			itemp=iresno(i)+1
		      else
			itemp=iresno(i)
		      endif
		        if(b(i,1).eq.' N ') itemp = iresno(na(i,l2))
			if(iresno(j).eq.itemp)then
			  mb(i,1)=i
			  mb(i,2)=j
c			  write(31,*)i,j
			endif
		    endif
		  enddo
c		  write(31,*) (mb(i,j1),j1=1,2), (b(i,j1),j1=1,2)
		  if(mb(i,2).ne.i.and.mb(i,2).ne.nb(i,1).and.
     &		    mb(i,2).ne.nb(i,2).and.mb(i,2).ne.nb(i,3).and.
     &		    mb(i,2).ne.nb(i,4).and.mb(i,2).ne.nb(i,5).and.
     &		    mb(i,2).ne.nb(i,6))then
			nb(i,kb(i))=mb(i,2)
			bn(i,kb(i))=b(i,2)
			kb(i)=kb(i)+1
		  endif
	       endif
	     enddo
	   endif
	 enddo
	enddo
c	write(31,*) i,(nb(i,j1),j1=1,kb(i)-1)
c	write(31,*) i,(bn(i,j1),j1=1,kb(i)-1)
C************************to find third neighbours***************************
	kc(i)=1
	do l3=1,kb(i)-1
	rewind 2
	    do i3=1,maxbnd
	      read(2,20)tcode,scode,iatom
		if(tcode.eq.rename(nb(i,l3)))then
	  	  do l2=1,iatom
		    read(2,30) (c(i,j),j=1,2)
		    if(bn(i,l3).eq.c(i,1))then
		      do j=1,natom
			if(atname(j).eq.c(i,2).and.c(i,2).ne.'    ')then
			  if(c(i,1).eq.' C  '.and.c(i,2).eq.' N  ')then
				itemp=iresno(i)+1
			  else
				itemp=iresno(i)
			  endif
		          if(c(i,1).eq.' N  ') itemp = iresno(nb(i,l3))
		          if(c(i,1).eq.' CA ') itemp = iresno(nb(i,l3))
			  if(iresno(j).eq.itemp)then
				mc(i,1)=i
				mc(i,2)=j
c				write(31,*)i,j
			  endif
			endif
		      enddo
c	write(31,*) (mc(i,j1),j1=1,2), (c(i,j1),j1=1,2)
c	if(mc(i,2).ne.i.and.mc(i,2).ne.na(i,1))then
	if(mc(i,2).ne.i.and.mc(i,2).ne.na(i,1).and.mc(i,2).ne.na(i,2)
     & .and.mc(i,2).ne.na(i,3).and.mc(i,2).ne.na(i,4).and.
     & mc(i,2).ne.na(i,5).and.mc(i,2).gt.i.and.mc(i,2).ne.nc(i,1).
     & and.mc(i,2).ne.nc(i,2).and.mc(i,2).ne.nc(i,3).and.mc(i,2).ne.
     & nc(i,4).and.mc(i,2).ne.nc(i,5))then
		nc(i,kc(i))=mc(i,2)
		cn(i,kc(i))=c(i,2)
		kc(i)=kc(i)+1
	endif
		endif
	      enddo
	    endif
	  enddo
	enddo
c	write(31,*) i,(nc(i,j1),j1=1,kc(i)-1)
c	write(31,*) i,(cn(i,j1),j1=1,kc(i)-1)
c****************************************************
      enddo
c
c	open(unit=4,file='mols.inp',status='unknown')
c	do i=1,natom
c	write(4,40)i,(ka(i)-1),(na(i,j1),j1=1,ka(i)-1)
c	write(31,*) i,(ka(i)-1),(na(i,j1),j1=1,ka(i)-1)
cc	write(31,*) i,(ka(i)-1),(an(i,j1),j1=1,ka(i)-1)
c	write(4,40)i,(kb(i)-1),(nb(i,j1),j1=1,kb(i)-1)
c	write(31,*) i,(kb(i)-1),(nb(i,j1),j1=1,kb(i)-1)
cc	write(31,*) i,(kb(i)-1),(bn(i,j1),j1=1,kb(i)-1)
c	write(4,40)i,(kc(i)-1),(nc(i,j1),j1=1,kc(i)-1)
c	write(31,*) i,(kc(i)-1),(nc(i,j1),j1=1,kc(i)-1)
cc	write(31,*) i,(kc(i)-1),(cn(i,j1),j1=1,kc(i)-1)
c	enddo
c	close(unit=4)
c
c
c	________________________________________________________________
c*********************to find fourth and above neighbours**********************
	do i = 1,natom
	  kd(i)=1
	  do j = i+1,natom
	    do j1=1,ka(i)-1
	      if(na(i,j1).eq.j) goto 111
	    enddo
	    do j2=1,kb(i)-1
	      if(nb(i,j2).eq.j) goto 111
	    enddo
	    do j3=1,kc(i)-1
	      if(nc(i,j3).eq.j) goto 111
	    enddo
 	    nd(i,kd(i))=j
c	    write(31,*) i,j
	    kd(i)=kd(i)+1
111	  enddo
c	  write(6,*) i,(nd(i,j1),j1=1,kd(i)-1)
	enddo
c*********************to count the atom pairs fourth & above neighbours********
	do i=1,natom
	  ke(i)=1
	  ne(i,ke(i))=nd(i,1)
	  do j=1,kd(i)-1
	    if(nd(i,j).ne.(nd(i,j-1)+1).and.j.gt.1)then
	      ke(i)=ke(i)+1
	      ne(i,ke(i))=nd(i,j-1)
	      ke(i)=ke(i)+1
	      ne(i,ke(i))=nd(i,j)
	    endif
	  enddo
	  ke(i)=ke(i)+1
	  ne(i,ke(i))=nd(i,kd(i)-1)
	  if(ne(i,1).eq.0) ke(i)=0
	enddo
c
c	do i=1,natom
c	write(31,*)i,ke(i)/2,(ne(i,j),j=1,ke(i))
c	enddo
c
c	open(unit=4,file='mols1.inp',status='unknown')
c	do i=1,natom
c	write(4,50) i,float(ke(i)/2),float(kc(i)-1)
c	write(4,60) (ne(i,j),j=1,ke(i)),(nc(i,j),j=1,kc(i)-1)
c	enddo
c	close(unit=4)
c
c***************************************************************
c	to exclude the atompairs in around peptide bond
c	and in the rings(1-5 pairs)
C**************************************************************
	do i=1,natom
	k15(i)=0
c	write(31,*) (ne(i,j),j=1,ke(i))
	do j=1,ke(i)
	j1=ne(i,j)
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HH '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' OH '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HH '.and.iresno(i).eq.iresno(j1)) goto 131

	if(rename(i).eq.'TYR '.and.atname(i).eq.' CD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HE1'.and.
     &iresno(i).eq.iresno(j1))  goto 151
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' CZ '.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' OH '.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'TYR '.and.atname(i).eq.' HH '.and.
     &iresno(i).eq.iresno(j1)) goto 151 

	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HZ '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ '.and.iresno(i).eq.iresno(j1)) goto 131

	if(rename(i).eq.'PHE '.and.atname(i).eq.' CD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HE1'.and.
     &iresno(i).eq.iresno(j1))  goto 151
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' CZ '.and.
     &iresno(i).eq.iresno(j1)) goto 151 
	if(rename(i).eq.'PHE '.and.atname(i).eq.' HZ '.and.
     &iresno(i).eq.iresno(j1)) goto 151 

	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE3'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HH2'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ2'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HZ3'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.' HH2'.and.iresno(i).eq.iresno(j1)) goto 131

	if(rename(i).eq.'TRP '.and.atname(i).eq.' CD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HD1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' NE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CE3'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HE3'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CZ2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HZ2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CZ3'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HZ3'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' CH2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'TRP '.and.atname(i).eq.' HH2'.and.
     &iresno(i).eq.iresno(j1)) goto 151

	if(rename(i).eq.'HIS '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE1'.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.' HE2'.and.iresno(i).eq.iresno(j1)) goto 131

	if(rename(i).eq.'HIS '.and.atname(i).eq.' HD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' ND1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HD2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' CE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HE1'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' NE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151
	if(rename(i).eq.'HIS '.and.atname(i).eq.' HE2'.and.
     &iresno(i).eq.iresno(j1)) goto 151

	if(rename(j1).eq.'PRO '.and.atname(i).eq.' C  '
     &.and.(iresno(i)+1).eq.iresno(j1)) goto 151
	if(rename(j1).eq.'PRO '.and.atname(i).eq.' O  '
     &.and.(iresno(i)+1).eq.iresno(j1)) goto 151
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CA '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CA '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'2HD '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HB '
     &.and.atname(j1).eq.'1HD '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'2HD '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'2HB '
     &.and.atname(j1).eq.'1HD '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CG '
     &.and.atname(j1).eq.'1HG '.and.iresno(i).eq.iresno(j1)) goto 131
	if(rename(i).eq.'PRO '.and.atname(i).eq.'1HG '
     &.and.atname(j1).eq.'2HG '.and.iresno(i).eq.iresno(j1)) goto 131
c
	goto 161
c
151	do j2=i+1,natom
	  if(atname(j2).eq.' C  ') then
		k15(i)=k15(i)+1
		n15(i,k15(i))=j2
		k15(i)=k15(i)+1
		n15(i,k15(i))=natom
		goto 141
	  endif
	enddo
c
161	k15(i)=k15(i)+1
	n15(i,k15(i))=ne(i,j)
c	
131	enddo
c	write(31,*) (n15(i,j),j=1,k15(i)),atname(i)
141	enddo
	
c******************************************************************
c	to exclude the atompairs in around peptide bond
c	and in the rings(in 1-4 pairs) and pickup needed 1-4 pairs
C******************************************************************
	do i = 1, natom
	k14(i)=1
	do l3=1,kb(i)-1
	rewind 2
	do i3=1,maxbnd
	read(2,20)tcode,scode,iatom
	if(tcode.eq.rename(nb(i,l3)))then
	do l2=1,iatom
	read(2,30) (dn(i,j),j=1,2)
	if(bn(i,l3).eq.dn(i,1))then
	do j=1,natom
c	if(atname(i).eq.' CA '.and.atname(j).eq.' CA '.and.
c     &(iresno(i)+1).eq.iresno(j)) goto 121
c	if(atname(i).eq.' CA '.and.atname(j).eq.' H  '.and.
c     &(iresno(i)+1).eq.iresno(j)) goto 121
c	if(atname(i).eq.' O  '.and.atname(j).eq.' CA '.and.
c     &(iresno(i)+1).eq.iresno(j)) goto 121
c	if(atname(i).eq.' O  '.and.atname(j).eq.' H  '.and.
c     &(iresno(i)+1).eq.iresno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE1'.and.atname(j).eq.' HE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE1'.and.atname(j).eq.' HH '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CE2'.and.atname(j).eq.' HH '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TYR ') goto 121
	if(atname(i).eq.' CD1'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'TYR '.and.
     &iresno(i).eq.iresno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'PHE ') goto 121

	if(atname(i).eq.' CD1'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CE1'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CE2'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CZ '.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HZ '.and.rename(i).eq.'PHE '.and.
     &iresno(i).eq.iresno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' HD1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' NE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' CE3'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CE2 '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE3'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' NE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CZ3'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HZ3'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'TRP ') goto 121

	if(atname(i).eq.' CD1'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD1'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' NE1'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CE2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CE3'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE3'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CZ2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HZ2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CZ3'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HZ3'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CH2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HH2'.and.rename(i).eq.'TRP '.and.
     &iresno(i).eq.iresno(j)) goto 121

	if(atname(i).eq.' CB '.and.atname(j).eq.' CE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' NE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CB '.and.atname(j).eq.' HD2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' NE2 '.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' HE2'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121
	if(atname(i).eq.' CG '.and.atname(j).eq.' CE1'.and.
     &iresno(i).eq.iresno(j).and.rename(i).eq.'HIS ') goto 121

	if(atname(i).eq.' ND1'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CE1'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE1'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' NE2'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HE2'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CD2'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' HD2'.and.rename(i).eq.'HIS '.and.
     &iresno(i).eq.iresno(j)) goto 121

	if(atname(i).eq.' C  '.and.rename(j).eq.'PRO ') goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'1HB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'2HB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.' CG '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.' CB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'1HG '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' N  '.and.
     &atname(j).eq.'2HG '.and.iresno(i).eq.iresno(j)) goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'1HB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'2HB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.' CG '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.' CB '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'1HG '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' HA '.and.
     &atname(j).eq.'2HG '.and.iresno(i).eq.iresno(j)) goto 121

	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.'1HD '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.'2HD '.and.iresno(i).eq.iresno(j)) goto 121
	if(rename(i).eq.'PRO '.and.atname(i).eq.' CB '.and.
     &atname(j).eq.' CD '.and.iresno(i).eq.iresno(j)) goto 121

	if(atname(i).eq.' CA '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'1HB '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'2HB '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CG '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'1HG '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'2HG '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.' CD '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'1HD '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
	if(atname(i).eq.'2HD '.and.rename(i).eq.'PRO '.and.
     &iresno(i).eq.iresno(j)) goto 121
c
	if(atname(j).eq.dn(i,2).and.dn(i,2).ne.'    ')then
	  if(dn(i,1).eq.' C  '.and.dn(i,2).eq.' N  ')then
	    itemp=iresno(i)+1
	  else
	    itemp=iresno(i)
	  endif
          if(dn(i,1).eq.' N  ') itemp = iresno(nb(i,l3))
          if(dn(i,1).eq.' CA ') itemp = iresno(nb(i,l3))
	  if(iresno(j).eq.itemp)then
	    md(i,1)=i
	    md(i,2)=j
c	    write(31,*)i,j
	  endif
	endif
121	enddo
c	write(31,*) (md(i,j1),j1=1,2), (dn(i,j1),j1=1,2)
c	if(md(i,2).ne.i.and.md(i,2).ne.na(i,1))then
	if(md(i,2).ne.i.and.md(i,2).ne.na(i,1).and.md(i,2).ne.na(i,2)
     &.and.md(i,2).ne.na(i,3).and.md(i,2).ne.na(i,4).and.
     &md(i,2).ne.na(i,5).and.md(i,2).gt.i.and.md(i,2).ne.n14(i,1).
     &and.md(i,2).ne.n14(i,2).and.md(i,2).ne.n14(i,3).and.
     &md(i,2).ne.n14(i,4).and.md(i,2).ne.n14(i,5))then
	n14(i,k14(i))=md(i,2)
	m14(i,k14(i))=dn(i,2)
	k14(i)=k14(i)+1
	endif
c
	endif
	enddo
	endif
	enddo
	enddo
c	write(31,*) i,(n14(i,j1),j1=1,k14(i)-1) 
c	write(31,*) i,(m14(i,j1),j1=1,k14(i)-1)
	enddo
	close(unit=2)
c
c***********************************************************************
c	modify here for ecepp.
	do i = 1,natom
	do j = i+1,natom
	if(rename(i).eq.'TYR'.and.iresno(i).eq.iresno(j).and.
     & atname(i).eq.' CA '.or.rename(i).eq.'TYR'.and.iresno(i)
     &.eq.iresno(j).and.atname(i).eq.'1HB'.or.rename(i).eq.
     &'TYR'.and.iresno(i).eq.iresno(j).and.atname(i).eq.'2HB'
     &) then
	if(atname(j).eq.' HD1'.or.atname(j).eq.' HD2'.or.
     & atname(j).eq.' HD1'.or.atname(j).eq.' HD2'.or.
     & atname(j).eq.' CE1'.or.atname(j).eq.' HE1'.or.
     & atname(j).eq.' CE2'.or.atname(j).eq.' HE2'.or.
     & atname(j).eq.' CZ '.or.atname(j).eq.' OH ') then
c	write (6,*) i,(n14(i,j1),j1=1,k14(i)-1)
  	n14(i,k14(i)) = j
	k14(i) = k14(i) + 1
c	print *,i,(n14(i,j1),j1=1,k14(i)-1)	
	endif
	else if(rename(i).eq.'PHE'.and.iresno(i).eq.iresno(j).and.
     & atname(i).eq.' CA '.or.rename(i).eq.'PHE'.and.iresno(i)
     &.eq.iresno(j).and.atname(i).eq.'1HB'.or.rename(i).eq.
     &'PHE'.and.iresno(i).eq.iresno(j).and.atname(i).eq.'2HB'
     &) then
	if(atname(j).eq.' HD1'.or.atname(j).eq.' HD2'.or.
     & atname(j).eq.' HD1'.or.atname(j).eq.' HD2'.or.
     & atname(j).eq.' CE1'.or.atname(j).eq.' HE1'.or.
     & atname(j).eq.' CE2'.or.atname(j).eq.' HE2'.or.
     & atname(j).eq.' CZ '.or.atname(j).eq.' HZ') then
c	write (6,*) i,(n14(i,j1),j1=1,k14(i)-1)
  	n14(i,k14(i)) = j
	k14(i) = k14(i) + 1
c	print *,i,(n14(i,j1),j1=1,k14(i)-1)	
	endif
	else if(rename(i).eq.'HIS'.and.iresno(i).eq.iresno(j).and.
     & atname(i).eq.' CA '.or.rename(i).eq.'HIS'.and.iresno(i)
     &.eq.iresno(j).and.atname(i).eq.'1HB'.or.rename(i).eq.
     &'HIS'.and.iresno(i).eq.iresno(j).and.atname(i).eq.'2HB'
     &) then
	if(atname(j).eq.' HD2'.or.atname(j).eq.' CE1'.or.atname(j)
     &.eq.' HE1'.or.atname(j).eq.' NE2'.or.atname(j).eq.' HE2') then
c	write (6,*) i,(n14(i,j1),j1=1,k14(i)-1)
  	n14(i,k14(i)) = j
	k14(i) = k14(i) + 1
c	print *,i,(n14(i,j1),j1=1,k14(i)-1)	
	endif
	else if(rename(i).eq.'TRP'.and.iresno(i).eq.iresno(j).and.
     & atname(i).eq.' CA '.or.rename(i).eq.'TRP'.and.iresno(i)
     &.eq.iresno(j).and.atname(i).eq.'1HB'.or.rename(i).eq.
     &'TRP'.and.iresno(i).eq.iresno(j).and.atname(i).eq.'2HB'
     &) then
	if(atname(j).eq.' HD1'.or.atname(j).eq.' NE1'.or.atname(j)
     &.eq.' HE1'.or.atname(j).eq.' CE2'.or.atname(j).eq.' CE3'.or.
     &atname(j).eq.' HE3'.or.atname(j).eq.' CZ2'.or.atname(j).eq.' HZ2'
     &.or.atname(j).eq.' CZ3'.or.atname(j).eq.' HZ3'.or.atname(j).eq.
     &' CH2'.or.atname(j).eq.' HH2') then
c	write (6,*) i,(n14(i,j1),j1=1,k14(i)-1)
  	n14(i,k14(i)) = j
	k14(i) = k14(i) + 1
c	print *,i,(n14(i,j1),j1=1,k14(i)-1)	
	endif
	else if(rename(i+5).eq.'PRO'.and.iresno(i+5).eq.iresno(j)
     &.and.atname(i).eq.' CA '.and.rename(i).ne.'PRO'
     &.or.rename(j).eq.'PRO'.and.iresno(i+1)
     &.eq.iresno(j).and.atname(i).eq.' O  ')then
cccc	print *,i,j,atname(i),atname(j),rename(i),rename(j)
	if(atname(j).eq.'1HB '.or.atname(j).eq.'2HB '.or.atname(j)
     &.eq.' CG '.or.atname(j).eq.'1HG '.or.atname(j).eq.'2HG '
     &.or.atname(j).eq.'1HD '.or.atname(j).eq.'2HD '
     &.or.atname(j).eq.' HA '.or.atname(j).eq.' C  '
     &.or.atname(j).eq.' CB ') then
cccc	print *,i,j,atname(i),atname(j),rename(i),rename(j)
c	write (6,*) i,(n14(i,j1),j1=1,k14(i)-1)
  	n14(i,k14(i)) = j
	k14(i) = k14(i) + 1
cccc	print *,i,(n14(i,j1),j1=1,k14(i)-1)	
	endif
	endif
  	enddo
	enddo
c
c******remove 1-4 pairs from the 1-5 pairs of rings*******************
c
c	do i = 83,83
c	print *,(n14(i,l),l=1,k14(i)-1)
cc	print *,(n15(i,k),k=1,k15(i))
c	enddo

	do i = 1,natom
	jj = 1
	jc = 1
	do j = 1,k15(i)/2
c	print *,j,jj,n15(i,jj),n15(i,jj+1)
	do k = n15(i,jj),n15(i,jj+1)
	do l = 1,k14(i)-1
	if (k.eq.n14(i,l)) go to 133
	enddo
	nt(i,jc) = k
c	print *,nt(i,jc)
	jc = jc + 1
133	enddo
	jj = jj + 2
	enddo
	kt(i) = jc
	enddo
c
	do i=1,natom
	  kf15(i)=1
	  nf15(i,kf15(i))=nt(i,1)
	  do j=1,kt(i)-1
	    if(nt(i,j).ne.(nt(i,j-1)+1).and.j.gt.1)then
	      kf15(i)=kf15(i)+1
	      nf15(i,kf15(i))=nt(i,j-1)
	      kf15(i)=kf15(i)+1
	      nf15(i,kf15(i))=nt(i,j)
	    endif
	  enddo
	  kf15(i)=kf15(i)+1
	  nf15(i,kf15(i))=nt(i,kt(i)-1)
	  if(nf15(i,1).eq.0) kf15(i)=0
	enddo
c     _________________________________________________________________________
c
c     WRITES 1-4,1-5 INTERACTIONS,NON-BONDED ENERGY PARAMETERS AND PARTIAL
c     CHARGE OF ATOMS, WHICH ARE APPLICABLE FOR ECEPP FORCE FIELD HERE
cc    _________________________________________________________________________
	open(unit=4,file=if1,status='unknown')
	open(unit=8,file='ATOMTYPE.lib',status='old')
	do i=1,natom
	if(i.eq.1)then
	 e1 = 0.0760
	if(rename(i).eq.'GLY') e1 = 0.080
	tt = 0.0
	aty = 4.0
	 write(4,100) i,e1,aty,tt,float(kf15(i)/2),float(k14(i)-1)
	 write(4,60) (nf15(i,j),j=1,kf15(i)),(n14(i,j),j=1,k14(i)-1)
	 go to 171
	endif
	rewind(8)
	do j=1,366
	  read(8,90) tcode,scode,natm
	  if(tcode.eq.rename(i))then
	    do j1=1,natm
	     read(8,80) at1,at2,at3,e2,e1,iaty
	if(atname(i).eq.' N  '.and.iresno(i).eq.1) e1 = -0.332
	if(atname(i).eq.' H  '.and.iresno(i).eq.1) e1 =  0.076
	if(atname(i).eq.' C  '.and.iresno(i).eq.lres) e1 = 0.517
	if(atname(i).eq.' O  '.and.iresno(i).eq.lres) e1 = -0.351
	if(tcode.eq.'GLY'.and.atname(i).eq.' N  '.and.
     &  iresno(i).eq.1) e1 = -0.328
	if(tcode.eq.'GLY'.and.atname(i).eq.' H  '.and.
     &  iresno(i).eq.1) e1 =  0.080
	if(tcode.eq.'GLY'.and.atname(i).eq.' C  '.and.
     &  iresno(i).eq.lres) e1 =  0.518
	if(tcode.eq.'GLY'.and.atname(i).eq.' O  '.and.
     &  iresno(i).eq.lres) e1 = -0.347
	if(rename(i-1).eq.'GLY'.and.atname(i).eq.' O1 '.and.
     &  iresno(i).eq.lres.and.tcode.eq.' OH') e1 = -0.335
	if(rename(i-2).eq.'GLY'.and.atname(i).eq.' HO '.and.
     &  iresno(i).eq.lres.and.tcode.eq.' OH') e1 =  0.230
	aty = float(iaty)
	     if(at1.eq.atname(i))then
	       if(iresno(i).eq.1.and.atname(i).eq.' N  ') at3 = 'NT'
	       write(4,100) i,e1,aty,tt,float(kf15(i)/2),float(k14(i)-1)
	       write(4,60) (nf15(i,j5),j5=1,kf15(i)),(n14(i,j4),
     &         j4=1,k14(i)-1)
	       go to 171
	     endif
	    enddo
	  endif
	enddo
 171	enddo
	close(unit=8)
ccc	close(unit=4)
c*******************write dihedral angles*************************************
	ntor = 0
	do i=1,lres
	rewind 10
	open(unit=10,file='DIHEDS.lib',status='old')
	do j=1,111
	read(10,90) tcode,scode,natm
	if(scode.eq.fseq(i))then
	if(iopt.eq.1) natm = 2
	do j1=1,natm
	id1= 0
	id2= 0
	id3= 0
	id4= 0
	read(10,110) d1,d2,d3,d4
	do j2=1,natom
	if(d1.eq.atname(j2).and.iresno(j2).eq.i) id1=j2
	if(d2.eq.atname(j2).and.iresno(j2).eq.i) id2=j2
	if(d3.eq.atname(j2).and.iresno(j2).eq.i) id3=j2
	if(j1.le.2) then
	  id4=natom
	else
	  if(d4.eq.atname(j2).and.iresno(j2).eq.i) id4=j2
	endif
	if(id1.ne.0.and.id2.ne.0.and.id3.ne.0.and.id4.ne.0)goto 333
	enddo
333	ntor = ntor + 1
	write(4,40) id1, id2, id3, id4
	write(31,*) id1, id2, id3, id4
	enddo
	endif
	enddo
	enddo
	write(31,*)'No. of Parameters(dihedral angles) = ',ntor
	write(*,*)'No. of Parameters(dihedral angles) = ',ntor

c******************write hydrogen bond parameters*****************************
	close(unit=4)
c*********************calculate aij & cij for non-bonded interaction***********
        conv= 332.0

        do i=1,mxtyat
          atpl(i)=atpl(i)/100
          efel(i)=efel(i)/100
          emin(i)=emin(i)/1000
          rmin(i)=rmin(i)/100
        enddo

        do i=1,mxtyat
          ri=rmin(i)
          ai=atpl(i)
          aei=sqrt(ai/efel(i))
cc          aic=ai/ehm             !!  ICM
cc          do j=i,mxtyat          !!  -"-
          aic=ai*ehm                 !!  comment for ICM:
          cij(i,i)=aic*ai/(aei+aei)  !!        -"-
          as=.5*cij(i,i)*ri**6
          aij(i,i)=as                 !!
          a14(i,i)=.5*as              !!
          do j=i+1,mxtyat            !!
            aj=atpl(j)
c _______ Constant for 6-12 attractive term (Slater-Kirkwood formula)
            cs=aic*aj/(aei+sqrt(aj/efel(j)))
            cij(i,j)=cs
            cij(j,i)=cs
c ____________________________ repulsive term (form. 3 & 6 of ref 2)
            rij=.5*(ri+rmin(j))
            as=.5*cs*rij**6
            aij(i,j)=as
            aij(j,i)=as
            as=.5*as
            a14(i,j)=as
            a14(j,i)=as
          enddo
c          do(i)=do_s(i)
c          ac(i)=ac_s(i)
        enddo

c +++++++++++++++++++++++++++++++++
        cij(1,1)=45.5d0
        aij(1,1)=14090.0d0
        cij(2,2)=45.5d0
        aij(2,2)=14380.0d0
        cij(3,3)=45.5d0
        aij(3,3)=8420.0d0
        cij(4,4)=45.5d0
        aij(4,4)=8420.0d0
        cij(5,5)=45.5d0
        aij(5,5)=11680.0d0
        cij(6,6)=45.5d0
        aij(6,6)=14380.0d0
        cij(7,7)=370.5
        aij(7,7)=906100.0
        cij(8,8)=766.6
        aij(8,8)=1049000.0
        cij(9,9)=509.5
        aij(9,9)=653600.0
        cij(10,10)=217.2
        aij(10,10)=125600.0
        cij(11,11)=369.0
        aij(11,11)=170200.0
        cij(12,12)=217.2
        aij(12,12)=125600.0
        cij(13,13)=401.3
        aij(13,13)=375200.0
        cij(14,14)=401.3
        aij(14,14)=375200.0
        cij(15,15)=401.3
        aij(15,15)=375200.0
        cij(16,16)=2274.4
        aij(16,16)=5809000.0
        cij(17,17)=45.5
        aij(17,17)=5340.0
        cij(18,18)=370.5
        aij(18,18)=909000.0
c +++++++++++++++++++++++++++++++++
        do i=1,mxtyat
          a14(i,i)=.5*aij(i,i)
        enddo
c +++++++++++++++++++++++++++++++++
c*****identify hydrogen bond pairs ******************************************************************************
c	hb-donors:3,4,5,6
c	hb-acceptors:10,11,12,15,16
        do i=1,mxtyat
          do(i)=do_f(i)
          ac(i)=ac_f(i)
        enddo
        do i=1,mxhbdo
          do j=1,mxhbac
            hbc(i,j)=chb_s(i,j)
            hba(i,j)=ahb_s(i,j)
          enddo
        enddo

      do i=1,mxtyat
        do j=1,mxtyat
          ihb(i,j)=0
          chb(i,j)=0.
          ahb(i,j)=0.
        enddo
      enddo

      iac=0
      ido=0
      do i=1,mxtyat
        if (do(i)) then
          ido=ido+1
          jac=0
          do j=1,i-1
            if (ac(j)) then
              jac=jac+1
              ihb(i,j)=1
              ihb(j,i)=-1
              chb(i,j)=hbc(ido,jac)
              chb(j,i)=chb(i,j)
              ahb(i,j)=hba(ido,jac)
              ahb(j,i)=ahb(i,j)
            endif
          enddo
        elseif (ac(i)) then
          iac=iac+1
          jdo=0
          do j=1,i-1
            if (do(j)) then
              jdo=jdo+1
              ihb(i,j)=-1
              ihb(j,i)=1
              chb(i,j)=hbc(jdo,iac)
              chb(j,i)=chb(i,j)
              ahb(i,j)=hba(jdo,iac)
              ahb(j,i)=ahb(i,j)
            endif
          enddo
        endif
      enddo
c

c******************************************************************************
c600	stop
600	return
	end

