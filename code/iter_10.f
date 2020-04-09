      subroutine iter(nt,Tq,TD,T3,dbeta,nodes,s,T0,ns,ne,nc,ni,nr)

      integer*8 nt, it, T0(nodes)
      integer s(nodes), Tq(nodes), T3(nodes), TD(nodes)
      double precision dbeta
      double precision dran_u, u
c
c 25-03: slope of deaths in IB: 0.241 ->>  0.252
c      
cc      print*,ns,ne,nc,ni,nr

      do it=nt+1, nt+nodes
        ix=i_dran(nodes)
	iy=i_dran(nodes)
	do while(ix.eq.iy)
	  iy=i_dran(nodes)
	enddo

c	
	if(s(ix)+s(iy).ne.0) then
	  if(s(ix).eq.0) then
	    if(s(iy).eq.2.or.s(iy).eq.3) then
	      u=dran_u()
	      if(u.lt.dbeta) then
	        s(ix)=1
		ns=ns-1
		ne=ne+1
		T0(ix)=it
	      endif
	    endif

	  elseif(s(ix).eq.1) then
           if(it-t0(ix).gt.Tq(ix)) then
	      nr=nr+1
	      ne=ne-1
	      s(ix)=4
            elseif(it-t0(ix).gt.TD(ix)) then
	      ni=ni+1
	      ne=ne-1
	      s(ix)=3
            elseif(it-t0(ix).gt.T3(ix)) then
	      nc=nc+1
	      ne=ne-1
	      s(ix)=2
c	      print*,ix,s(ix), T0(ix), it, T3(ix), Td(3)
c	      stop
 	    endif  

	  elseif(s(ix).eq.2) then
           if(it-t0(ix).gt.Tq(ix)) then
	      nr=nr+1
	      nc=nc-1
	      s(ix)=4
            elseif(it-t0(ix).gt.TD(ix)) then
	      ni=ni+1
	      nc=nc-1
	      s(ix)=3
 	    endif  

	  elseif(s(ix).eq.3) then
            if(it-t0(ix).gt.Tq(ix)) then
	      nr=nr+1
	      ni=ni-1
	      s(ix)=4
	    endif  
	  endif 


	endif

c       print*,it,ns,ne,nc,ni,nr
      enddo
c      stop
      nt=it-1
      

      end
