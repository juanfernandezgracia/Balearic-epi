c Total pop: 1095426
c municipios = 67
      parameter(nodes=1095426,np=67,nrun=100)
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
     
      integer Tq, T30
      double precision dbeta, dmu
      
      Tq=2*14*nodes
      T30=2*30*nodes
      dbeta=0.3d0
      dmu=.1d0
c---------------------------
      
      do i=1, np
        npobl1(i)=0
        npobl2(i)=0
	do j=1, np
	  nl(i,j)=0
	enddo	  
      enddo
      
c      open(1,file='comm105.csv')
      open(1,file='mobility.csv')
      read(1,*)
      
      do i=1, 200
        read(1,*,END=100) n1, n2, n12
	npobl1(n1+1)=npobl1(n1+1) + n12
	npobl2(n2+1)=npobl2(n2+1) + n12
	nl(n1+1,n2+1)=n12
c        read(1,*,END=100) chr1, n11
c	npobl1(i)=npobl1(i) + n11
c	npobl2(i)=npobl2(i) + n11
c	nl(i,i)=n11
      enddo
100   continue
c      print*,npobl1
c      print*,npobl2
       ntmp1=0
       ntmp2=0
c       do i=1, np
c         ntmp1=ntmp1+ npobl1(i)
c         ntmp2=ntmp2 + npobl2(i)
c	 print*,i, npobl1(i),npobl2(i), ntmp1, ntmp2
c       enddo
cc      stop
      
      index1(0)=0
      index1(1)=0
      index2(0)=0
      index2(1)=0
      do i=2, np
        index1(i)=index1(i-1) + npobl1(i-1)
        index2(i)=index2(i-1) + npobl2(i-1)
      enddo
c      print*, index2(36), index1(0)
      
c      print*,index1
c      print*,index2
c      stop
      
      ix=0      
      do i=1, np
        do j=1, np
	  do k=1, nl(i,j)
	    ix=ix+1
	    index2(j)=index2(j)+1
	    l2(ix)=index2(j)
c	    print*,ix,l2(ix)
	  enddo
	enddo
      enddo
c      print*,8698, l2(8698)
cc      
      ix=0      
      do j=1, np
        do i=1, np
	  do k=1, nl(i,j)
	    ix=ix+1
	    index1(i)=index1(i)+1
	    l1(ix)=index1(i)
c	    if(j.eq.36) then
c	      print*,'c',i,j,ix,l1(ix), index1(i)
c	      stop
c	    endif  
	  enddo
	enddo
c	print*,j,ix
      enddo
c      print*, 'a',index2(35), index1(1), index2(1), l1(386526)
      
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
	dbeta=.1d0*dble(ib)*.1d0
	
	do i=1, 1000
          hist_r(i)=0
          hist_i(i)=0
          hist_r2(i)=0.d0
          hist_i2(i)=0.d0
	enddo	
 
 	dbeta=0.1d0
	nr_a=0
	do ir=1, nrun
	  do i=1, nodes
	    s1(i)=0
	    s2(i)=0
	    T1(i)=0
	    T2(i)=0
	  enddo

	  ix=i_dran(nodes)
	  s1(ix)=1
	  s2(l2(ix))=s1(ix)
	  ns=nodes-1
	  ne=1
	  ni=0
	  nr=0
c---------
	  it=0
	  do while(ne+ni.ne.0)
c	  do irep=1, 10
	    do i=1, np
              call iter(it,Tq,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
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
c	  print*,it, dbeta,ntmp1, ntmp2, ns,ne,ni,nr
c	  stop

	    do i=1, np
              call iter(it,Tq,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr,npobl2(i)
	    enddo
c	  print*,it
	    do i=1, nodes
c 	  print*,'a',it,l1(i), l2(i)
             s1(l1(i))=s2(i)
c	  print*,'b',it
              T1(l1(i))=T2(i)
	    enddo
c	  write(10,*) it/2/nodes, ns,ne,ni,nr

	    ind=it/2/nodes
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
	    enddo
	    write(8,*) it/2/nodes, dbeta, (nr_tmp(i), i=1, np)
 	  enddo
c---------
	    do i=1, np
              nr_tmp(i)=0
 	      do j= index1(i-1)+ 1, index1(i)
		if(s1(j).eq.3) nr_tmp(i) = nr_tmp(i) + 1
	      enddo
	    enddo
	    write(11,*) it/2/nodes, dbeta, (nr_tmp(i), i=1, np)
c---	    
	  nr_a=nr_a + nr
	  ind=it/2/nodes
	  print*,ind, dbeta,ns,ne,ni,nr

c	  do i=1, np
c            nr_tmp=0
c 	    do j= index1(i-1)+ 1, index1(i)
c	      if(s1(j).eq.3) nr_tmp = nr_tmp + 1
c	    enddo
c	    write(8,*) it/2/nodes, dbeta, i, nr_tmp
c	  enddo
          do i=ind+1, 1000
	    hist_r(i)=hist_r(i) + nr
	    hist_r2(i)=hist_r2(i) + dble(nr)*dble(nr)
	  enddo
	  
	enddo
	
	print*,dbeta, dble(nr_a)/dble(nrun)
	do i=1, 1000
	  av_r= dble(hist_r(i))/dble(nrun)
	  av2_r=dble(hist_r2(i))/dble(nrun)
	  sd_r=dsqrt(av2_r -av_r*av_r)
	 
	 av_i= dble(hist_i(i))/dble(nrun)
	  av2_i=dble(hist_i2(i))/dble(nrun)
	  sd_i=dsqrt(av2_i -av_i*av_i)
	  
c	  if(hist_r(i).ne.0) write(9,*) i, dble(hist_r(i))/dble(nrun)
c     &  , sd_r, dble(hist_i(i))/dble(nrun), sd_i
	  if(hist_r(i).ne.0) write(9,*) i, av_r, sd_r, av_i, sd_i
	enddo
      enddo	  
      
      end
