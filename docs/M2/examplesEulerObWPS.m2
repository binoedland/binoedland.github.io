loadPackage ("EulerObstructionWPS",Reload=>true)

--gives correct answer of 69
dualDegreeSurface(2,3,5)

--gives correct answers 2688,40
dualDegreeThreefold(2,3,5)
dualDegreeThreefold(6,10,15)


--Generate examples in article appendix: isolated singularities
l = {};
for k from 1 to 10 do (
    for m from k to 10 do (
	for n from m to 10 do(
	    if gcd(k,m) == 1 and gcd(k,n)==1 and gcd(m,n)==1 then (
		e1=euObThreefold(k,m,n);
		e2=euObThreefold(m,n,k);
	    	e3=euObThreefold(n,k,m);
		r1=RSV(k,m,n);
	    	r2=RSV(m,n,k);
		r3=RSV(n,k,m);
		l = l | {{k,m,n,e1,e2,e3,r1,r2,r3}};
	       	);
	    );
	);
    );
l = matrix l

--Generate examples in article appendix: Non-isolated singularities
l = {};
for k from 1 to 6 do (
    for m from k to 6 do (
	for n from m to 6 do(
	    if ((gcd(k,m) == 1 and gcd(k,n)==1 and gcd(m,n)==1) == false) and (gcd(k,m,n) ==1) then (
		e1=euObThreefold(k,m,n);
		e2=euObThreefold(m,n,k);
	    	e3=euObThreefold(n,k,m);
		r1=RSV(k,m,n);
	    	r2=RSV(m,n,k);
		r3=RSV(n,k,m);
		l = l | {{k,m,n,e1,e2,e3,r1,r2,r3}};
	       	);
	    );
	);
    );
l = matrix l
