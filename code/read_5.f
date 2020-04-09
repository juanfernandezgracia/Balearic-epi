c Total pop: 1095426
c municipios = 67
      parameter(nodes=1095426,np=67,nrun=4)
c s2(l2(i)) = s1(i)
c s1(l1(i)) = s2(i)
c
      integer*8 it, T1(nodes), T2(nodes)
      integer s1(nodes), s2(nodes)
      integer l1(nodes), l2(nodes)
      integer npobl1(np), npobl2(np), nl(np,np)
      integer index1(0:np), index2(0:np)
      character*50 chr1
      
      integer nr_a, nr_tmp(np)
      integer hist_tinf(1000)
      integer hist_r(1000), hist_i(1000)
      double precision hist_r2(1000), hist_i2(1000)
      double precision av_r, av2_r, av_i, av2_i, sd_r, sd_i
     
      integer Tq_1(nodes), T3_1(nodes), Toff
      integer Tq_2(nodes), T3_2(nodes)
      integer Tinc
      double precision dbeta, dbeta1, dbeta2, dbeta3, dbeta4, dbeta5
      
      character*22 infile1, infile2, infile3, infile4, infile5
      CHARACTER*10 buffer
 
      call dran_ini(918988765)

      CALL getarg(1, buffer)
      read(buffer,*) Tinc 
      
      do i=1, nodes
c       T3_1(i)=2*(Tinc+i_dran(2)-1)*nodes
        T3_1(i)=2*Tinc*nodes
      enddo
      
      CALL getarg(2, buffer)
      read(buffer,*) Tdis 

      do i=1, nodes
c        Tq_1(i)=2*(Tdis+Tinc+i_dran(2)-1)*nodes
        Tq_1(i)=2*(Tdis+Tinc)*nodes
      enddo

      CALL getarg(3, buffer)
      read(buffer,*) dbeta1 

      CALL getarg(4, buffer)
      read(buffer,*) dbeta2

      CALL getarg(5, buffer)
      read(buffer,*) dbeta3

      CALL getarg(6, buffer)
      read(buffer,*) dbeta4

      dbeta1=dbeta1
      dbeta2=dbeta2
      dbeta3=dbeta3
      dbeta4=dbeta4
      dbeta5=0.d0
      
      
      Toff=20
      iphase1=Toff + 15
      iphase2=Toff + 22
      iphase3=Toff + 30
      iphase4=Toff + 57
      
      infile1='../SEIR/mobility.csv'
      infile2='../SEIR/mobility.csv'
      infile3='../SEIR/mobility.csv'
      infile4='../SEIR/mobility.csv'
      infile5='../SEIR/mobility.csv'
      
c      write(8,'(a)') "#irun CODIGOINE beta time recover" 
      write(9,'(a)') "#time average_recover sd_R average_infected sd_I" 
c      write(11,'(a)') "#time beta nr_final_per_city" 
      write(12,'(a)') "#time susceptible exposed infected recover" 
      write(13,'(a)') "#time number_cases_inititated_at_time" 
      
c------------------------

      do ib=1, 1
c	dbeta=.1d0*dble(ib)*.1d0
	
	do i=1, 1000
          hist_r(i)=0
          hist_i(i)=0
          hist_r2(i)=0.d0
          hist_i2(i)=0.d0
	  hist_tinf(i)=0
	enddo	
 
	nr_a=0
	do ir=1, nrun
c---Phase I
	  dbeta=dbeta1
 
          call meta(infile1,nodes,np,npobl1,npobl2,index1,index2,l1,l2)

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i)
	      stop
	    endif  
	  enddo

c--------------
	  do i=1, nodes
	    s1(i)=0
	    s2(i)=0
	    T1(i)=0
	    T2(i)=0
	    Tq_2(l2(i))=Tq_1(i)
	    T3_2(l2(i))=T3_1(i)
	  enddo

cPalma
	  ix=index1(42) + i_dran(npobl1(43))
c	  ix=i_dran(nodes)
	  s1(ix)=1
	  s2(l2(ix))=s1(ix)
	  ns=nodes-1
	  ne=1
	  ni=0
	  nr=0
	  it=0
c---
	 do nit=1, iphase1
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq_1(index1(i-1) +1),
     &T3_1(index1(i-1)+1),dbeta,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo
c check!!!
            ntmp1=0
            ntmp2=0
	    do i=1, nodes
 		if(s1(i).eq.0) ntmp1 = ntmp1 + 1
 		if(s2(i).eq.0) ntmp2 = ntmp2 + 1
	    enddo

	    do i=1, np
              call iter(it,Tq_2(index2(i-1) +1),
     &T3_2(index2(i-1)+1),dbeta,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
              s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind-Toff, ns, ne, ni, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ni
	    hist_i2(ind)=hist_i2(ind) + dble(ni)*dble(ni)
c--c	    print*,ind, ni, hist_i(ind)

c: Daily stats
c
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
c	    enddo
c	    do i=1, np
c-----c	      write(8,*) ir, i, dbeta, it/2/nodes,  nr_tmp(i)
	    enddo
 	  enddo
c---Phase II
c---------------------------   
          dbeta=dbeta2   

          call meta(infile2,nodes,np,npobl1,npobl2,index1,index2,l1,l2)

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i)
	      stop
	    endif  
	  enddo

c------------------------
	  do nit=iphase1+1, iphase2
c	  do irep=1, 10
	    do i=1, np
               call iter(it,Tq_1(index1(i-1) +1),
     &T3_1(index1(i-1)+1),dbeta,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
cc             call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
cc     &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo
c check!!!
            ntmp1=0
            ntmp2=0
	    do i=1, nodes
 		if(s1(i).eq.0) ntmp1 = ntmp1 + 1
 		if(s2(i).eq.0) ntmp2 = ntmp2 + 1
	    enddo

	    do i=1, np
              call iter(it,Tq_2(index2(i-1) +1),
     &T3_2(index2(i-1)+1),dbeta,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
cc              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
cc     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo
            
	    
	    ind=it/2/nodes
	    write(12,*) ind-Toff, ns, ne, ni, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ni)*dble(ni)
c--c	    print*,ind, ni, hist_i(ind)
c: Daily stats
c
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
c	    enddo
c	    do i=1, np
c----c	      write(8,*) ir, i, dbeta, it/2/nodes,  nr_tmp(i)
	    enddo
 	  enddo

c---Phase III
c---------------------------   
          dbeta=dbeta3   

          call meta(infile3,nodes,np,npobl1,npobl2,index1,index2,l1,l2)

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i)
	      stop
	    endif  
	  enddo

c------------------------
	  do nit=iphase2+1, iphase3
c	  do irep=1, 10
	    do i=1, np
               call iter(it,Tq_1(index1(i-1) +1),
     &T3_1(index1(i-1)+1),dbeta,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
cc             call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
cc     &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo
c check!!!
            ntmp1=0
            ntmp2=0
	    do i=1, nodes
 		if(s1(i).eq.0) ntmp1 = ntmp1 + 1
 		if(s2(i).eq.0) ntmp2 = ntmp2 + 1
	    enddo

	    do i=1, np
              call iter(it,Tq_2(index2(i-1) +1),
     &T3_2(index2(i-1)+1),dbeta,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
cc              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
cc     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo
            
	    
	    ind=it/2/nodes
	    write(12,*) ind-Toff, ns, ne, ni, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ni)*dble(ni)
c--c	    print*,ind, ni, hist_i(ind)
c: Daily stats
c
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
c	    enddo
c	    do i=1, np
c----c	      write(8,*) ir, i, dbeta, it/2/nodes,  nr_tmp(i)
	    enddo
 	  enddo

c---Phase IV
c---
	  dbeta=dbeta4

          call meta(infile4,nodes,np,npobl1,npobl2,index1,index2,l1,l2)

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i)
	      stop
	    endif  
	  enddo

c------------------------
         do nit=iphase3+1, iphase4
c	  do while(ne+ni.ne.0)
c	  do irep=1, 10
	    do i=1, np
               call iter(it,Tq_1(index1(i-1) +1),
     &T3_1(index1(i-1)+1),dbeta,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
cc             call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
cc     &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo
c check!!!
            ntmp1=0
            ntmp2=0
	    do i=1, nodes
 		if(s1(i).eq.0) ntmp1 = ntmp1 + 1
 		if(s2(i).eq.0) ntmp2 = ntmp2 + 1
	    enddo

	    do i=1, np
              call iter(it,Tq_2(index2(i-1) +1),
     &T3_2(index2(i-1)+1),dbeta,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
cc              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
cc     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind-Toff, ns, ne, ni, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ni
	    hist_i2(ind)=hist_i2(ind) + dble(ni)*dble(ni)
c: Daily stats
c
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
c	    enddo
c	    do i=1, np
c-----c	      write(8,*) ir, i, dbeta, it/2/nodes,  nr_tmp(i)
	    enddo
 	  enddo
c	  stop
cc	  goto 1000
c---------
c---Phase V
c---
	  dbeta=dbeta5

          call meta(infile5,nodes,np,npobl1,npobl2,index1,index2,l1,l2)

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i)
	      stop
	    endif  
	  enddo

c------------------------
	  do while(ne+ni.ne.0)
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq_1(index1(i-1) +1),
     &T3_1(index1(i-1)+1),dbeta,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
cc              call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
cc    &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo
c check!!!
            ntmp1=0
            ntmp2=0
	    do i=1, nodes
 		if(s1(i).eq.0) ntmp1 = ntmp1 + 1
 		if(s2(i).eq.0) ntmp2 = ntmp2 + 1
	    enddo

	    do i=1, np
              call iter(it,Tq_2(index2(i-1) +1),
     &T3_2(index2(i-1)+1),dbeta,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
cc              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
cc     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind-Toff, ns, ne, ni, nr
	    
	    if(ind.le.1000) then 
	      hist_r(ind)=hist_r(ind) + nr
	      hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	      hist_i(ind)=hist_i(ind) + ni
	      hist_i2(ind)=hist_i2(ind) + dble(ni)*dble(ni)
	    endif
c: Daily stats
c
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
c	    enddo
c	    do i=1, np
c-----c	      write(8,*) ir, i, dbeta, it/2/nodes,  nr_tmp(i)
	    enddo
 	  enddo
c---------
c---------
c---------
1000	  do i=1, np
            nr_tmp(i)=0
 	    do j= index1(i-1)+ 1, index1(i)
	      if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	    enddo
	  enddo
c-----c	  write(11,*) it/2/nodes, dbeta, (nr_tmp(i), i=1, np)
c---	    
	  nr_a=nr_a + nr
	  ind=it/2/nodes
c....c	  print*,ind, dbeta,ns,ne,ni,nr

c          if(ind.gt.1000) then
c	    print*,'increase hist!!!'
c	    stop
c	  endif          
          if(ind.lt.1000) then
	    do i=ind+1, 1000
	      hist_r(i)=hist_r(i) + nr
	      hist_r2(i)=hist_r2(i) + dble(nr)*dble(nr)
	    enddo
	  endif
	  
	  do i=1, nodes
	    j=T1(i)/2/nodes
	    if(j.gt.0) hist_tinf(j) = hist_tinf(j) + 1
	  enddo
	enddo
	

	do i=1, 1000
	  av_r= dble(hist_r(i))/dble(nrun)
	  av2_r=dble(hist_r2(i))/dble(nrun)
	  sd_r=dsqrt(av2_r -av_r*av_r)

	  	 av_i= dble(hist_i(i))/dble(nrun)
	  av2_i=dble(hist_i2(i))/dble(nrun)
	  sd_i=dsqrt(av2_i -av_i*av_i)
	  
	  write(9,*) i-Toff, av_r, sd_r, av_i, sd_i
	  if(hist_tinf(i).ne.0) 
     &                 write(13,*) i-Toff, dble(hist_tinf(i))/dble(nrun)
	enddo
      enddo	  
      
      end
