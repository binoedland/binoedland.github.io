newPackage("EulerObstructionWPS",
    Version => "0.1",
    Date => "January 17, 2017",
    Authors => {{Name=> "Bernt Ivar Utstoel Noedland"}},
    PackageExports => {"Polyhedra"},
    Headline => "Computes local euler obstructions and dual degrees of weighted projective surfaces and 3-folds"
    );

export{"RSV",
       "euObThreefold",
       "euObSurface",
       "HJFraction",
       "dualDegreeSurface",
       "dualDegreeThreefold"
    }

dualDegreeSurface = method();
dualDegreeSurface (ZZ,ZZ,ZZ) := (k,m,n) -> (
    a := 0;
    b := 0;
    c := 0;
    for i from 1 to k do (
	if (m+i*n)%k == 0 then (
	    a=i;
	    break;
	    );
	); 
    for i from 1 to m do (
	if (n+i*k)%m == 0 then (
	    b=i;
	    break;
	    );
	); 
    for i from 1 to n do (
	if (k+i*m)%n == 0 then (
	    c=i;
	    break;
	    );
	);
    deg := 3*k*m*n-2*(k+m+n)+euObSurface(k,k-a)+euObSurface(m,m-b)+euObSurface(n,n-c);
    deg);


dualDegreeThreefold = method();
dualDegreeThreefold (ZZ,ZZ,ZZ) := (k,m,n) -> (
    delta := lcm(k,m,n);
    vol := sub(delta^3/(k*m*n),ZZ);
    area := sub((delta^2*(k+m+n+1))/(k*m*n),ZZ);
    xylength :=sub((delta*gcd(k,m))/(k*m),ZZ);
    xzlength :=sub((delta*gcd(k,n))/(k*n),ZZ);
    yzlength :=sub((delta*gcd(m,n))/(m*n),ZZ);
    l := sub((delta*(k*m+k*n+m*n))/(k*m*n),ZZ)+xylength*euObSurface(gcd(k,m),n%(gcd(k,m)))+xzlength*euObSurface(gcd(k,n),m%(gcd(k,n)))+yzlength*euObSurface(gcd(m,n),k%(gcd(m,n)));
    vert := 1 + euObThreefold(k,m,n)+euObThreefold(m,n,k)+euObThreefold(n,k,m);
    deg := 4*vol-3*area+2*l-vert;
    deg);

RSV = method();
RSV (ZZ,ZZ,ZZ) := (k,m,n) -> (
    delta := lcm(k,m,n);
    A := transpose (matrix{{-1,0,0}});
    for x from 0 to (sub (delta/k -1,ZZ)) do (
	for y from 0 to (sub (delta/m,ZZ)) do (
	    for z from 0 to (sub (delta/n,ZZ)) do (
		if k*x+m*y+n*z <= delta then (
		    if k*x+(m-gcd(k,m))*y+(n-gcd(k,n))*z >= delta-k then (
			A = A | (matrix{{x},{y},{z}});
			);
		    );
		);
	    );
	);
    Q := convexHull A;
    B := A | (matrix{{sub(delta/k,ZZ)},{0},{0}});
    P := convexHull B;
    rsv := 6*(volume(P)-volume(Q));		     
    rsv);

euObThreefold = method();
euObThreefold (ZZ,ZZ,ZZ) := (k,m,n) -> (
    c := 0;
    l := sub(k/(gcd(k,m)*gcd(k,n)),ZZ);
    for i from 1 to l do (
	if (m*gcd(k,n)+i*n)%k == 0 then (
	    c=i;
	    break;
	    );
	);
    eu := RSV(k,m,n);
    eu = eu - 6 +euObSurface(sub((k/gcd(k,n)),ZZ),sub(n/gcd(k,n),ZZ)%sub(k/gcd(k,n),ZZ))+euObSurface(sub(k/gcd(k,m),ZZ),sub(m/gcd(k,m),ZZ)%sub(k/gcd(k,m),ZZ))+euObSurface(sub(k/(gcd(k,m)*gcd(k,n)),ZZ),sub((k-c*gcd(n,k))/(gcd(k,m)*gcd(k,n)),ZZ));
    eu = eu +1 +euObSurface(gcd(k,n),m%gcd(n,k))+euObSurface(gcd(k,m),n%gcd(k,m));  
    sub(eu,ZZ));


euObSurface = method();
euObSurface (ZZ,ZZ) := (d,k) -> (
    HJ := HJFraction(d,k);
    eu := 0;
    for i from 0 to (length HJ -1) do (
	eu = eu +2 - HJ_i;
	);
    eu);

HJFraction = method();
HJFraction (ZZ,ZZ) := (d,k) -> (
    HJ := {};
    if d == 1 then (HJ = {1})
    else (
    	while k != 0 do (
	    q :=  ceiling(d/k);
	    HJ = HJ | {q};
	    r := q*k-d;
	    d = k;
	    k =r;
	    );
	);
    HJ);






