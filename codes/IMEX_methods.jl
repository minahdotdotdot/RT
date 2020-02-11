include("MMT.jl")

#ARK3(2)4L[2]SA /KC(3,4,3)
Ae = zeros(3,3); Ae[1,1]=1767732205903/2027836641118;
Ae[2,1:2]=[5535828885825/10492691773637;
 788022342437/10882634858940];
Ae[3,1:3]=[6485989280629/16251701735622;
-4246266847089/9704473918619;
10755448449292/10357097424841];

d=1767732205903/4055673282236;

Ai = zeros(3,3); Ai[1,1] = d;
Ai[2,1:2] = [2746238789719/10658868560708;
-640167445237/6845629431997];
Ai[3,1:3] = [1471266399579/7840856788654;
-4482444167858/7529755066697;
11266239266428/11593286722821];

b = [1471266399579/7840856788654;
-4482444167858/7529755066697;
11266239266428/11593286722821;
1767732205903/4055673282236];

c = [0; 1767732205903/2027836641118; 3/5; 1];

ARK3 = IMEXTableau(Ae, Ai, b, c, d);

#ARK4(3)6L[2]SA /KC(5,6,4)
Ae = zeros(5,5); Ae[1,1]=1/2;
Ae[2,1:2] = [13861/62500; 6889/62500];
Ae[3,1:3] = [-116923316275/2393684061468;
-2731218467317/15368042101831;
9408046702089/11113171139209];
Ae[4,1:4] = [-451086348788/2902428689909;
-2682348792572/7519795681897;
12662868775082/11960479115383;
3355817975965/11060851509271];
Ae[5,1:5] = [647845179188/3216320057751;
73281519250/8382639484533;
552539513391/3454668386233;
3354512671639/8306763924573;
4040/17871];

d = 1/4;

Ai = zeros(5,5); Ai[1,1]=d;
Ai[2,1:2] = [8611/62500; -1743/31250];
Ai[3,1:3] = [5012029/34652500;
-654441/2922500;
174375/388108];
Ai[4,1:4] = [15267082809/155376265600;
-71443401/120774400;
730878875/902184768;
2285395/8070912];
Ai[5,1:5] = [82889/524892;
0;
15625/83664;
69875/102672;
-2260/8211];

b = [82889/524892; 
0; 
15625/83664; 
69875/102672;
-2260/8211;
1/4]

c = [1/2; 83/250; 31/50; 17/20; 1];
ARK4 = IMEXTableau(Ae, Ai, b, c, d);

IMEXdict= Dict( 
	"ARK3" => ARK3, 
	"ARK4" => ARK4
	);

