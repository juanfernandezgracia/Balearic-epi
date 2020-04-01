c Total pop: 1095426
c municipios = 67
      parameter(nodes=1095426,np=67,nrun=1)
c s2(l2(i)) = s1(i)
c s1(l1(i)) = s2(i)
c
      integer s1(nodes), s2(nodes), T1(nodes), T2(nodes)
      integer l1(nodes), l2(nodes)
      integer npobl1(np), npobl2(np), nl(np,np)
      integer index1(0:np), index2(0:np)
      character*50 chr1
      
      integer nr_a, nr_tmp(np)
      integer hist_r(1000), hist_i(1000)
      double precision hist_r2(1000), hist_i2(1000)
      double precision av_r, av2_r, av_i, av2_i, sd_r, sd_i
     
      integer Tq, T3
      integer Tinc
      double precision dbeta, dmu
      
      CHARACTER*10 buffer
 
      CALL getarg(1, buffer)
      read(buffer,*) Tinc 
      T3=2*Tinc*nodes
      
      CALL getarg(2, buffer)
      read(buffer,*) Tdis 
      Tq=2*dis*nodes

      CALL getarg(3, buffer)
      read(buffer,*) dbeta 

c      dbeta=0.29d0
c      dmu=.1d0
      
      iphase1=20 + 8
      iphase2=20 + 30
      iphase3=20 + 45
      
c      write(8,*) "#irun CODIGOINE beta time recover" 
      write(9,*) "#time average_recover sd_R average_infected sd_I" 
c      write(11,*) "#time beta nr_final_per_city" 
      write(12,*) "#time susceptible infected recover" 
      
c---------------------------      
      do i=1, np
        npobl1(i)=0
        npobl2(i)=0
	do j=1, np
	  nl(i,j)=0
	enddo	  
      enddo
      
      open(1,file='mobility.csv')
      read(1,*)
      
      do i=1, 200
        read(1,*,END=100) n1, n2, n12
	npobl1(n1+1)=npobl1(n1+1) + n12
	npobl2(n2+1)=npobl2(n2+1) + n12
	nl(n1+1,n2+1)=n12
      enddo
100   continue
      close(1)
      
      index1(0)=0
      index1(1)=0
      index2(0)=0
      index2(1)=0
      do i=2, np
        index1(i)=index1(i-1) + npobl1(i-1)
        index2(i)=index2(i-1) + npobl2(i-1)
      enddo
      
      ix=0      
      do i=1, np
        do j=1, np
	  do k=1, nl(i,j)
	    ix=ix+1
	    index2(j)=index2(j)+1
	    l2(ix)=index2(j)
	  enddo
	enddo
      enddo
cc      
      ix=0      
      do j=1, np
        do i=1, np
	  do k=1, nl(i,j)
	    ix=ix+1
	    index1(i)=index1(i)+1
	    l1(ix)=index1(i)
	  enddo
	enddo
      enddo
      
      do i=1, nodes
        if(l1(l2(i)).ne.i) then
	  print*,'s',i,l1(l2(i)), l2(i), nl(1, 36)
	  stop
	endif  
      enddo
c      stop
c------------------------
      call dran_ini(918988765)

      do ib=1, 1
c	dbeta=.1d0*dble(ib)*.1d0
	
	do i=1, 1000
          hist_r(i)=0
          hist_i(i)=0
          hist_r2(i)=0.d0
          hist_i2(i)=0.d0
	enddo	
 
c 	dbeta=0.1d0
	nr_a=0
	do ir=1, nrun
	  do i=1, nodes
	    s1(i)=0
	    s2(i)=0
	    T1(i)=0
	    T2(i)=0
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
c---Phase I
	  do nit=1, iphase1
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
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
              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
              s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind, ns, ne+ne, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ne +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ne +ni)*dble(ne +ni)

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
c          dbeta=0.15d0   
	  do i=1, np
            npobl1(i)=0
            npobl2(i)=0
	    do j=1, np
	      nl(i,j)=0
	    enddo	  
	  enddo

	  open(1,file='mobility50.csv')
	  read(1,*)

	  do i=1, 200
            read(1,*,END=200) n1, n2, n12
	    npobl1(n1+1)=npobl1(n1+1) + n12
	    npobl2(n2+1)=npobl2(n2+1) + n12
	    nl(n1+1,n2+1)=n12
	  enddo
200       continue
	  close(1)

	  index1(0)=0
	  index1(1)=0
	  index2(0)=0
	  index2(1)=0
	  do i=2, np
            index1(i)=index1(i-1) + npobl1(i-1)
            index2(i)=index2(i-1) + npobl2(i-1)
	  enddo

	  ix=0      
	  do i=1, np
            do j=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index2(j)=index2(j)+1
		l2(ix)=index2(j)
	      enddo
	    enddo
	  enddo
cc      
	  ix=0      
	  do j=1, np
            do i=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index1(i)=index1(i)+1
		l1(ix)=index1(i)
	      enddo
	    enddo
	  enddo

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i), nl(1, 36)
	      stop
	    endif  
	  enddo

c------------------------
	  do nit=iphase1+1, iphase2
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
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
              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo
            
	    
	    ind=it/2/nodes
	    write(12,*) ind, ns, ne+ne, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ne +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ne +ni)*dble(ne +ni)
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
c---
c	  dbeta=0.1d0
	  do i=1, np
            npobl1(i)=0
            npobl2(i)=0
	    do j=1, np
	      nl(i,j)=0
	    enddo	  
	  enddo

	  open(1,file='mobility30.csv')
	  read(1,*)

	  do i=1, 200
            read(1,*,END=300) n1, n2, n12
	    npobl1(n1+1)=npobl1(n1+1) + n12
	    npobl2(n2+1)=npobl2(n2+1) + n12
	    nl(n1+1,n2+1)=n12
	  enddo
300       continue
	  close(1)

	  index1(0)=0
	  index1(1)=0
	  index2(0)=0
	  index2(1)=0
	  do i=2, np
            index1(i)=index1(i-1) + npobl1(i-1)
            index2(i)=index2(i-1) + npobl2(i-1)
	  enddo

	  ix=0      
	  do i=1, np
            do j=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index2(j)=index2(j)+1
		l2(ix)=index2(j)
	      enddo
	    enddo
	  enddo
cc      
	  ix=0      
	  do j=1, np
            do i=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index1(i)=index1(i)+1
		l1(ix)=index1(i)
	      enddo
	    enddo
	  enddo

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i), nl(1, 36)
	      stop
	    endif  
	  enddo

c------------------------
          do nit=iphase2+1, iphase3
c	  do while(ne+ni.ne.0)
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
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
              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind, ns, ne+ne, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ne +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ne +ni)*dble(ne +ni)
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
cc	  goto 1000
c---------
c---Phase IV
c---
	  do i=1, np
            npobl1(i)=0
            npobl2(i)=0
	    do j=1, np
	      nl(i,j)=0
	    enddo	  
	  enddo

	  open(1,file='mobility.csv')
	  read(1,*)

	  do i=1, 200
            read(1,*,END=400) n1, n2, n12
	    npobl1(n1+1)=npobl1(n1+1) + n12
	    npobl2(n2+1)=npobl2(n2+1) + n12
	    nl(n1+1,n2+1)=n12
	  enddo
400       continue
	  close(1)

	  index1(0)=0
	  index1(1)=0
	  index2(0)=0
	  index2(1)=0
	  do i=2, np
            index1(i)=index1(i-1) + npobl1(i-1)
            index2(i)=index2(i-1) + npobl2(i-1)
	  enddo

	  ix=0      
	  do i=1, np
            do j=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index2(j)=index2(j)+1
		l2(ix)=index2(j)
	      enddo
	    enddo
	  enddo
cc      
	  ix=0      
	  do j=1, np
            do i=1, np
	      do k=1, nl(i,j)
		ix=ix+1
		index1(i)=index1(i)+1
		l1(ix)=index1(i)
	      enddo
	    enddo
	  enddo

	  do i=1, nodes
            if(l1(l2(i)).ne.i) then
	      print*,'s',i,l1(l2(i)), l2(i), nl(1, 36)
	      stop
	    endif  
	  enddo

c------------------------
	  do while(ne+ni.ne.0)
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq,T3,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
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
              call iter(it,Tq,T3,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
	    enddo
	    do i=1, nodes
             s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo

	    ind=it/2/nodes
	    write(12,*) ind, ns, ne+ne, nr
	    
	    hist_r(ind)=hist_r(ind) + nr
	    hist_r2(ind)=hist_r2(ind) + dble(nr)*dble(nr)
	    hist_i(ind)=hist_i(ind) + ne +ni
	    hist_i2(ind)=hist_i2(ind) + dble(ne +ni)*dble(ne +ni)
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
	  print*,ind, dbeta,ns,ne,ni,nr

          if(ind.gt.1000) then
	    print*,'increase hist!!!'
	    stop
	  endif          
          do i=ind+1, 300
	    hist_r(i)=hist_r(i) + nr
	    hist_r2(i)=hist_r2(i) + dble(nr)*dble(nr)
	  enddo
	  
	enddo
	

	do i=1, 1000
	  av_r= dble(hist_r(i))/dble(nrun)
	  av2_r=dble(hist_r2(i))/dble(nrun)
	  sd_r=dsqrt(av2_r -av_r*av_r)

	  	 av_i= dble(hist_i(i))/dble(nrun)
	  av2_i=dble(hist_i2(i))/dble(nrun)
	  sd_i=dsqrt(av2_i -av_i*av_i)
	  
	  if(hist_r(i).ne.0) write(9,*) i, av_r, sd_r, av_i, sd_i
	enddo
      enddo	  
      
      end
