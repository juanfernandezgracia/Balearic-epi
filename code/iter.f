      subroutine iter(nt,Tq,T3,dbeta,dmu,nodes,s,T0,ns,ne,ni,nr)

      integer s(nodes), T0(nodes), Tq, T3
      double precision dbeta, dmu
      double precision dran_u, u
c
c 25-03: slope of deaths in IB: 0.241 ->>  0.252
c      
      

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
	      if(it-t0(iy).gt.t3) then
		u=dran_u()
		if(u.lt.dbeta) then
	          s(ix)=1
		  ns=ns-1
		  ne=ne+1
		  T0(ix)=it
		endif
	      endif
	    endif

	  elseif(s(ix).eq.1) then
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
            if(it-t0(ix).gt.Tq) then
	      nr=nr+1
	      ni=ni-1
	      s(ix)=3
	    endif  
	  endif 


	endif

      enddo
      
      nt=it-1
      

      end
