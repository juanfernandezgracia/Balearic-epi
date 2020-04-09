      parameter(ndays=37)
      double precision dinf(ndays), di_acu(ndays)
      double precision a1,a 2, a3
      double precision dalfa_i, x2_i
      double precision dalfa_r, x2_r      
      double precision dm_i(-19:1000), dm_ir(-19:1000)
      
      open(1,file='fort.9')
      open(2,file='../SEIR/data_inf.dat')
      
      read(1,*)
      do i=1, 1000
        read(1,*,END=100) nday, a1, a2, a3
	if(a1+a3.ne.0.d0) then
          dm_ir(nday)=dlog(a1+a3)
 	else
	  dm_ir(nday)=-10.d0
	endif
	if(a3.ne.0.d0) then
	  dm_i(nday)=dlog(a3)
 	else
	  dm_i(nday)=-10.d0
	endif
      enddo
100   continue

      ndayi =1
      ndayf = ndays     
      do i=ndayi, ndayf
        read(2,*,END=200) n1, n2, n3, n4
	dinf(i)=dlog(dble(n2-n3-n4))
	di_acu(i)=dlog(dble(n2))
      enddo
200   continue
      close(2)

c--- Inf activos
c      ndayi_adj= 15
      ndayi_adj= 22
      d1=0.d0
      d2=0.d0
      do i=ndayi_adj, ndayf
        d1=d1 + dm_i(i)
	d2=d2 + dinf(i)
      enddo
      dalfa_i= (d2-d1)/dble(ndayf-ndayi_adj+1)
c      print*,'inf= ',dalfa_i
            
      do i=ndayi, ndayf
        write(21,*) i, dexp(dalfa_i+dm_i(i)), dexp(dinf(i))
      enddo
      
      d1=0.d0
      do i=ndayi, ndayf
        d1=d1 + (dalfa_i + dm_i(i)- dinf(i))**2.d0
      enddo
      x2_i=d1
c      print*,'inf Xi2= ', x2_i
      
c-----      

c---Inf acumulados
      ndayi_adj= 15
      d1=0.d0
      d2=0.d0
      do i=ndayi_adj, ndayf
        d1=d1 + dm_ir(i)
	d2=d2 + di_acu(i)
      enddo
      dalfa_r= (d2-d1)/dble(ndayf - ndayi_adj +1)
c      print*,'muertos= ',dalfa_r
            
      do i=ndayi_adj, ndayf
        write(20,*) i, dexp(dalfa_r+dm_ir(i)), dexp(di_acu(i))
      enddo
      
      d1=0.d0
      do i=ndayi_adj, ndayf
        d1=d1 + (dalfa_r + dm_ir(i)- di_acu(i))**2
c	print*,d1, dalfa_r*dm_ir(i), death(i)
      enddo
      x2_r=d1
c      print*,'mue Xi2= ', x2_r
c-----      
      
      print*, dexp(dalfa_r), x2_r, dexp(dalfa_i), X2_i
       
      end      
