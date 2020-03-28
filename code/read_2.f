      parameter(nodes=100000,np=2,nrun=10)
c s2(l2(i)) = s1(i)
c s1(l1(i)) = s2(i)
c
      integer s1(nodes), s2(nodes), T1(nodes), T2(nodes)
      integer l1(nodes), l2(nodes)
      integer npobl1(np), npobl2(np), nl(np,np)
      integer index1(0:np), index2(0:np)
      
      integer nr_a
      
      integer Tq, T30
      double precision dbeta, dmu
      
      Tq=2*14*nodes
      T30=2*30*nodes
      dbeta=0.5d0
      dmu=.1d0
c---------------------------
      
      do i=1, np
        npobl1(i)=0
        npobl2(i)=0
      enddo
      
      open(1,file='comm105.csv')
      
      do i=1, 100
        read(1,*,END=100) n1, n2, n12
	npobl1(n1)=npobl1(n1) + n12
	npobl2(n2)=npobl2(n2) + n12
	nl(n1,n2)=n12
      enddo
100   continue
cc      print*,npobl1
cc      print*,npobl2
      
      index1(0)=0
      index1(1)=0
      index2(0)=0
      index2(1)=0
      do i=2, np
        index1(i)=index1(i-1) + npobl1(i-1)
        index2(i)=index2(i-1) + npobl2(i-1)
      enddo
      
cc      print*,index1
cc      print*,index2
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
cc      
      ix=0      
      do i=1, np
        do j=1, np
	  do k=1, nl(j,i)
	    ix=ix+1
	    index1(j)=index1(j)+1
	    l1(ix)=index1(j)
c	    print*,ix,l1(ix)
	  enddo
	enddo
      enddo
      
c      do i=1, nodes
c        print*,i,l1(l2(i))
c      enddo
c------------------------
      call dran_ini(918988765)

      do ib=1, 10
	dbeta=.1d0*dble(ib)*.1d0
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

	  it=0
	  do while(ne+ni.ne.0)
	    do i=1, np
              call iter(it,Tq,dbeta,dmu,npobl1(i),s1(index1(i-1) +1),
     &T1(index1(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s2(l2(i))=s1(i)
              T2(l2(i))=T1(i)
	    enddo

	    do i=1, np
              call iter(it,Tq,dbeta,dmu,npobl2(i),s2(index2(i-1) +1),
     &T2(index2(i-1) +1),ns,ne,ni,nr)
c	  print*,it,ns,ne,ni,nr
	    enddo
	    do i=1, nodes
              s1(l1(i))=s2(i)
              T1(l1(i))=T2(i)
	    enddo
	  enddo
	  nr_a=nr_a + nr
c	  print*,it/nodes,dbeta,ns,ne,ni,nr, nr_a
	enddo
	print*,dbeta, dble(nr_a)/dble(nrun)
      enddo	  
      
      end
