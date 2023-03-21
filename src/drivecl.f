c     library name : drive.f   

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
c     License along with this library; if not, write to the Free Software
c     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
c
c
c     contact : n_gautham@hotmail.com  

	Program MAIN
	
	include 'mols.par'
	parameter(mxtyat = 18)
	common /par/ natom, ntor, nhb, ns, lres
	common /tor/ u0(maxpar),sn(maxpar),tn(maxpar),phi(maxpar)
	common /ecpp/ aij(mxtyat,mxtyat),cij(mxtyat,mxtyat),
     &  a14(mxtyat,mxtyat),ihb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),
     &  chb(mxtyat,mxtyat)
	common /ctl/iopt,iff,icint,fseq(maxres)
c
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
cr      WRITE (31,1)
cr      WRITE (*,1)
cr1     FORMAT (//
cr     *10X,'#########################################################'/
cr     *10X,'#                                                       #'/
cr    *10X,'#     T  H  E    M  O  L  S    P  R  O  G  R  A  M      #'/
cr     *10X,'#                    Version 1.0                        #'/
cr     *10X,'#                                                       #'/
cr     *10X,'#     Department of Crystallography and Biophysics,     #'/ 
cr     *10X,'#    University of Madras, Chennai - 600 025, INDIA     #'/
cr     *10X,'#                                                       #'/
cr     *10X,'#                                                       #'/
cr     *10X,'#########################################################'/)
	CALL DATE(DAY)
	CALL TIME(HOUR)
	
	
cr      WRITE (*,'(2(/25X,2A/))') 'Job started on ', DAY,
cr     &                               '      at ', HOUR
cr      WRITE (31,'(2(/25X,2A/))') 'Job started on ', DAY, 
cr     &                               '      at ', HOUR
c--------read inputs-------------------------------
cr	print *,'molecule name :  ',mname
	read(30,'(A)') seq
	read(30,*) ns
	read(30,*) iopt
	read(30,*) irun
	read(30,*) iff
	read(30,*) ico
	read(30,*) icnt
	read(30,*) rk1,kr2
c

c******************************************************
2	format(/a55/)
cr	write(31,2) '***********GENERATING INITIAL COORDINATES**************'
cr	write(*,2) '***********GENERATING INITIAL COORDINATES**************'
	call pdbgen(seq,fseq)
cr	print *,'pdbgen o.k'
c
c	print *,iopt,irun,iff,icnt
c	if(irun.eq.1) go to 919
c	write(31,2) '*************GENERATING MOLS PARAMETERS****************'
c	write(*,2) '*************GENERATING MOLS PARAMETERS****************'
c	if(iff.eq.1) call amppar(iopt)
c	if(iff.eq.2) call ecppar(iopt,fseq,itcount)
c*******Flushing standard output************
c	call flush(6)
c*******End of Flushing ********************
c	print *,'pargen o.k'
c*******Flushing standard output************
c	call flush(6)
c*******End of Flushing ********************
	call varinit
c	call conformation(0)

c	write(31,2) '********GENERATING MOLS OPTIMAL STRUCTURES*************'
c	write(*,2) '********GENERATING MOLS OPTIMAL STRUCTURES*************'
cr	call mols(fseq,iopt,iff)
c	print *,'mols o.k'
c*******Flushing standard output************
c	stop
	call flush(6)
c*******End of Flushing ********************
c	call conformation(1)
c	if(irun.eq.2) go to 919
c
c	write(31,2) '*******MINIMISING MOLS OPTIMAL STURUCTRES BY CG********'
c	write(*,2) '*******MINIMISING MOLS OPTIMAL STURUCTRES BY CG********'
c	call minimiz
c	print *,'minimiz o.k'
c*******Flushing standard output************
c	call flush(6)
c*******End of Flushing ********************
c	call conformation(2)
c	if(irun.eq.3) go to 919
c
	write(31,2) '**********CLUSTERING MINIMIZED STRUCTURES**************'
	write(*,2) '**********CLUSTERING MINIMIZED STRUCTURES**************'
	if(ico.eq.1) call cluster(icnt)
	if(ico.eq.2) call mclust(icnt,rk1,kr2)
c
	print *,'cluster o.k'
c
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************
c******************************************************
cr	write(31,2) '******************RUN INFORMATION*********************'
cr	write (31,*)'Output directory  :  ',path
cr	write (31,*)'Molecule name       :  ',mname
cr	write (31,*)'Sequence            :  ',seq
cr	write (31,*)'Number of Cycles    :  ',ns
c	If(iopt.eq.1) then
c	  write (31,*)'Search option       :  ','Back bone only'
c	 elseif(iopt.eq.2) then
c	  write (31,*)'Search option       :  ','Back bone + side chain'
c	 else
c	  write (31,*)'Search option       :  ','Rotomer option'
c	endif
c	if(irun.eq.1) then
c	  write (31,*)'Run option          :  ','Upto Pdbgen'
c	 elseif(irun.eq.2) then
c	  write (31,*)'Run option          :  ','Upto MOLS'
c	 elseif(irun.eq.3) then
c	  write (31,*)'Run option          :  ','Upto Minimiz'
c	 else
c	  write (31,*)'Run option          :  ','Upto Clustering'
c	endif
cr	if(iff.eq.1) then
cr	  write (31,*)'Force field         :  ','AMBER'
cr	 else
cr	  write (31,*)'Force field         :  ','ECEPP/3'
cr	endif
cr	if(ico.eq.1) then
cr	  write (31,*)'Cluster algorithm   :  ','K-Means algorithm'
cr	 else
cr	  write (31,*)'Cluster algorithm   :  ','Hiearchical algorithm'
cr	endif
cr	write (31,*)'Cluster interval    :  ',icnt
cr	if(ico.eq.2) then
cr	write (31,*)'Cluster cutoff      :  ',rk1
cr	If(kr2.eq.1) then
cr	  write (31,*)'Atoms for clustering:  ','All atoms'
cr	 elseif(kr2.eq.2) then
cr	  write (31,*)'Atoms for clustering:  ','Backbone atoms'
cr	 else
cr	  write (31,*)'Atoms for clustering:  ','Calfa atoms'
cr	endif
cr	endif
919      CALL DATE(DAY)
      CALL TIME(HOUR)
      WRITE (31,'(2(/25X,2A/))') 'Job ended on ', DAY, 
     &                               '      at ', HOUR
      WRITE (*,'(2(/25X,2A/))') 'Job ended on ', DAY,
     &                               '      at ', HOUR
	write (*,*)'cpu time = : ',etime(tt)
	write(31,*)'cpu time = : ',etime(tt)
c*******Flushing standard output************
	call flush(6)
c*******End of Flushing ********************

	stop
	end
