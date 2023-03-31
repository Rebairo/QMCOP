c     library name : drive.f   

c     Copyright (c) 2023      Abraham Rebairo, Sam Paul D and Stephen A

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
c     License along with this library; if not, write to the Free Software
c     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
c
c
c     contact : abrahamrebairo@hotmail.com  
c       
	Program MAIN
	
	include 'mols.par'
	parameter(mxtyat = 18)
	common /par/ natom, ntor, nhb, ns, lres
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
	common /ctl/iopt,iff,icint,fseq(maxres)
	common /sopt/isopt
	common /fnam/ if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
	common/mout/opmo(maxstr,maxpar),opmi(maxstr,maxpar),
     &              emo(maxstr),emi(maxstr)
	character seq*(maxres),mname*128,path*128
	character*128 if1,pf1,pf2,pf3,pf4,of1,of2,of3,of4,of0
c
      CHARACTER*8 HOUR
      CHARACTER*9 DAY
	real tt(2)
	open(unit=30,file='usermols.inp',status='old')
	read(30,'(A)') path
c	print *,'path name :  ',path
	read(30,*) mname
c	print *,'molecule name :  ',mname
c
c***********Out put generate file names**********************
	of0 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'inf.log'
	if1 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'inp.par'
	pf1 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'inp.pdb'
	pf2 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'mols.pdb'
	pf3 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'mini.pdb'
	pf4 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'cent.pdb'
	of1 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'mols.out'
	of2 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'mini.out'
	of3 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'clut.pdb'
	of4 = path(:LNBLNK(path)) // '/' // 
     &        mname(:LNBLNK(mname)) // '_' // 'cplt.pdb'
c	print *,if1,of0,of1,pf1,pf2,pf3,pf4,of1,of2,of3,of4
c************************************************************
	open(unit=31,file=of0,status='unknown')
c      WRITE (31,1)
c      WRITE (*,1)
c1     FORMAT (//
c     *10X,'#########################################################'/
c     *10X,'#                                                       #'/
c     *10X,'#     T  H  E    Q U A S I P O L Y   P R O G R A M      #'/
c     *10X,'#                    Version 1.0                        #'/
c     *10X,'#                                                       #'/
c     *10X,'#         Department of Nuclear Physics,           ,     #'/ 
c     *10X,'#    University of Madras, Chennai - 600 025, INDIA     #'/
c     *10X,'#                                                       #'/
c     *10X,'#                                                       #'/
c     *10X,'#########################################################'/)
	CALL DATE(DAY)
	CALL TIME(HOUR)
	
	
c      WRITE (*,'(2(/25X,2A/))') 'Job started on ', DAY,
c     &                               '      at ', HOUR
c      WRITE (31,'(2(/25X,2A/))') 'Job started on ', DAY, 
c     &                               '      at ', HOUR
c--------read inputs-------------------------------
!	print *,'molecule name :  ',mname
	read(30,'(A)') seq
!	read(30,*) ns
	read(30,*) iopt
	read(30,*) irun,isopt
	read(30,*) iff
	read(30,*) ico
	read(30,*) icnt
	read(30,*) rk1,kr2
c
        open(unit=134,file='mols_q.out',status='old')
                i=0
135             continue
                read(134,*,END=136)
                        i=i+1
                        
                goto 135
136             continue
                ns = i
       close(unit=134)
!        print *,ns

c******************************************************
2	format(/a55/)
	call pdbgen(seq,fseq)
	
        if(iff.eq.1) call amppar(iopt)
cr	if(iff.eq.2) call ecppar(iopt,fseq,itcount)
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
!	print *,'pargen o.k'
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
	call varinit
c	call conformation(0)

	call mols(fseq,iopt,iff)
cr	print *,'mols o.k'
c*******Flushing standard output************
c	stop
	call flush(6)
c*******End of Flushing ********************
	call conformation(1)
c*******Flushing standard output************
	call flush(6)
	write(31,2) '******************RUN INFORMATION*********************'
	write (31,*)'Output directory  :  ',path
	write (31,*)'Molecule name       :  ',mname
	write (31,*)'Sequence            :  ',seq
	write (31,*)'Number of Cycles    :  ',ns
	If(iopt.eq.1) then
	  write (31,*)'Search option       :  ','Back bone only'
	 elseif(iopt.eq.2) then
	  write (31,*)'Search option       :  ','Back bone + side chain'
	 else
	  write (31,*)'Search option       :  ','Rotomer option'
	endif
cr	if(irun.eq.1) then
	  write (31,*)'Run option          :  ','Upto Pdbgen'
cr	 elseif(irun.eq.2) then
	  write (31,*)'Run option          :  ','Upto MOLS'
cr	 elseif(irun.eq.3) then
	  write (31,*)'Run option          :  ','Upto Minimiz'
cr	 else
	  write (31,*)'Run option          :  ','Upto Clustering'
cr	endif
	if(iff.eq.1) then
	  write (31,*)'Force field         :  ','AMBER'
	 else
	  write (31,*)'Force field         :  ','ECEPP/3'
	endif
	if(ico.eq.1) then
	  write (31,*)'Cluster algorithm   :  ','K-Means algorithm'
	 else
	  write (31,*)'Cluster algorithm   :  ','Hiearchical algorithm'
	endif
	write (31,*)'Cluster interval    :  ',icnt
	if(ico.eq.2) then
	write (31,*)'Cluster cutoff      :  ',rk1
	If(kr2.eq.1) then
	  write (31,*)'Atoms for clustering:  ','All atoms'
	 elseif(kr2.eq.2) then
	  write (31,*)'Atoms for clustering:  ','Backbone atoms'
	 else
	  write (31,*)'Atoms for clustering:  ','Calfa atoms'
	endif
	endif
919      CALL DATE(DAY)
      CALL TIME(HOUR)
      WRITE (31,'(2(/25X,2A/))') 'Job ended on ', DAY, 
     &                               '      at ', HOUR
	write(31,*)'cpu time = : ',etime(tt)
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************

	stop
	end
