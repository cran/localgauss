../src/obfun_dv_dv.f90 : ../src/obfun.f90 ../src/obfun_dv.f90 
	tapenade -o obfun_dv -d -multi -root "lgobfun_dv" -vars "pp" -outvars "lld" ../src/obfun_dv.f90  -html
	sed 's/USE DIFFSIZES/ /' < obfun_dv_dv.f90 > tmp.f90
	sed 's/nbdirsmax/5/' < tmp.f90 > tmp1.f90
	sed 's/nbdirsmax/5/' < tmp1.f90 > tmp2.f90
	sed 's/nbdirsmax/5/' < tmp2.f90 > tmp3.f90
	sed 's/nbdirsmax/5/' < tmp3.f90 > tmp4.f90		
	sed 's/nbdirsmax/5/' < tmp4.f90 > obfun_dv_dv.f90
	cp obfun_dv_dv.f90 ../src/
	rm tmp*



../src/obfun_dv.f90 : ../src/obfun.f90
	tapenade -o obfun -d -multi -root "lgobfun" -vars "pp" -outvars "ll" ../src/obfun.f90 -html
	sed 's/REAL,/REAL(8),/' < obfun_dv.f90 > tmp.f90
	sed 's/REAL ::/REAL(8) ::/' < tmp.f90 > obfun_dv.f90
	sed 's/USE DIFFSIZES/ /' < obfun_dv.f90 > tmp.f90
	sed 's/nbdirsmax/5/' < tmp.f90 > tmp1.f90
	sed 's/nbdirsmax/5/' < tmp1.f90 > tmp2.f90
	sed 's/nbdirsmax/5/' < tmp2.f90 > tmp3.f90
	sed 's/nbdirsmax/5/' < tmp3.f90 > tmp4.f90		
	sed 's/nbdirsmax/5/' < tmp4.f90 > obfun_dv.f90
	cp obfun_dv.f90 ../src/


