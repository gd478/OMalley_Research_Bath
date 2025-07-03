PROGRAM PLOT
	!Alex Porter
	REAL :: time(22000), a(22000, 56)
	CHARACTER*50 :: nam, EU
	INTEGER :: nstep(22000), nument(22000), i, j, x
100 FORMAT (1i10 , 1e14.6, 1i10)
200 FORMAT (5e14.6)
300 FORMAT (1a50)
	
	!open and read input file
	OPEN (UNIT=1,FILE = 'STATIS', status = 'old')
	READ(1,300) nam
	PRINT*, nam
	READ(1,300) EU	
	PRINT*, EU
	DO i=1,22000
		READ(1,100) nstep(i), time(i), nument(i)
		!PRINT*, nstep(i), time(i), nument(i)
			READ(1,200) (a(i,j),j=1, nument(i))
			!PRINT*, a
	ENDDO
		
	!Write output file 
	OPEN(UNIT=2, FILE = 'TEMP.csv', STATUS = 'unknown')
	WRITE(2,*) ' nstep',',','step',',', 'time',',', 'TEMP',',', ' engcns',',', ' engsrc',',', ' engcpe'
	DO x = 1, 22000
		WRITE(2,*) x,',', nstep(x),',', time(x),',', a(x, 2),',', a(x, 1),',', a(x, 4),',', a(x, 5),','	
	ENDDO	

END PROGRAM PLOT