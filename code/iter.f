      subroutine iter(nt,Tq,dbeta,dmu,nodes,s,T0,ns,ne,ni,nr)

      integer s(nodes), T0(nodes), Tq, T30
      double precision dbeta, dmu
      double precision dran_u, u
c
c 25-03: slope of deaths in IB: 0.241 ->>  0.252
c      
      

cc      Tq=14*nodes
cc      T30=30*nodes
cc     dbeta=0.4d0
cc      dmu=.1d0
cc      dgamma=0.01d0
      
             
cc      do ib=1, 1
cc	dbeta=0.252d0 + 0.1d0*dble(ib)*.003d0
        
cc	dbeta=0.4d0
	
cc	nr_a=0
cc	nr_t=0
cc	n30=0

	do it=nt+1, nt+nodes
          ix=i_dran(nodes)
	  iy=i_dran(nodes)
	  do while(ix.eq.iy)
	    iy=i_dran(nodes)
	  enddo
c	
	  if(s(ix)+s(iy).ne.0) then
	    if(s(ix).eq.0) then
	      if(s(iy).eq.1.or.s(iy).eq.2) then
		u=dran_u()
		if(u.lt.dbeta) then
	          s(ix)=1
		  ns=ns-1
		  ne=ne+1
		  T0(ix)=it
		endif
	      endif

	    elseif(s(ix).eq.1) then
c		if(dran_u().lt.dgamma) then
              if(it-t0(ix).gt.Tq) then
		nr=nr+1
		ne=ne-1
		s(ix)=3
	      elseif(dran_u().lt.dmu) then 
		ni=ni+1
		ne=ne-1
		s(ix)=2	      
	      endif  

	    elseif(s(ix).eq.2) then
c		if(dran_u().lt.dgamma) then
              if(it-t0(ix).gt.Tq) then
		nr=nr+1
		ni=ni-1
		s(ix)=3
	      endif  
	    endif 
cc	if(ix.eq.1) write(8,*) it, ns, ne, ni, nr, ix, s(ix), iy, s(iy)


	  endif

cc	  if(mod(it,nodes).eq.0.and.nr.ne.0) then
cc           write(8,*) it/nodes,nr
cc	    write(10,*) it/nodes,ne+ni
cc         endif  
cc	  if(it.eq.t30) n30=n30 + nr 
	enddo
c      print*,dbeta,nr
cc	nr_a=nr_a+ nr
cc	nr_t=nr_t+ it

cc      enddo
      
      nt=it-1
      

      end
