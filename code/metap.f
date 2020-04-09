      subroutine meta(infile,nodes,np,npobl1,npobl2,index1,index2,l1,l2)
      
      integer l1(nodes), l2(nodes)
      integer npobl1(np), npobl2(np), nl(np,np)
      integer index1(0:np), index2(0:np)
      
            character*22 infile
	  
      do i=1, np
	npobl1(i)=0
	npobl2(i)=0
	do j=1, np
	  nl(i,j)=0
	enddo	  
      enddo

c      open(1,file='mobility.csv')
      open(1,file=infile)
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
	  print*,'s',i,l1(l2(i)), l2(i)
	  stop
	endif  
      enddo

      end

