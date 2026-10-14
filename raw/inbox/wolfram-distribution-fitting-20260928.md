# Wolfram computations for distribution fitting, 2026-09-28

Queries use synthetic distributions and equations, with no private observations sent.
The first natural-language moment query returned only its input; the reformulated
query below solved it. The exact-law formula is derived in the diagnostic report;
Wolfram independently checks its Jacobian, moments and numerical integrals.

## moment_equations

<result query='solve a*n=m, 2*a^2*n=v for a,n' url='https://www.wolframalpha.com/input?i=solve+a%2An%3Dm%2C+2%2Aa%5E2%2An%3Dv+for+a%2Cn'>

# Input interpretation

solve a n = m
2 a^2 n = v for a, n

# Results

a = 0 and v = 0 and m = 0

a!=0 and n = 0 and v = 0 and m = 0

a = v/(2 m) and n = (2 m^2)/v and m!=0 and v!=0

</result>

## symbolic_moments

Out[1]= {{{a -> v/(2*m), nu -> (2*m^2)/v}}, (b + Cos[b]*Sin[b])/(2*b), (12*b + 8*Sin[2*b] + Sin[4*b])/(32*b), 0.321750554396642193401404614358661319`18., 0.20483276469913345164919784755051894808`18.}

## exact_law_checks

Out[1]= {1, {1.93239932581490158454485996508943216472`18., 3.74134960950164143631942758241886627744`18., 0.96805809221753306340042937395554277971`18., 1.99616050043892477048314009221062436674`18.}, {{1, 0.5957907655987166618116296634457010110280643365579715618916`35.}, {4, 0.1262007995358370369957356399944611546344979403514466117562`35.}, {7, 0.0267936646320212407733599555596530869808753991209262142891`35.}, {10, 0.0057011906706149227420095412515882693385329264396309786131`35.}, {16, 0.0002597559459337312712960254323670991176789280823457781696`35.}}, {{1, 0.5168539955343700451718766814640245434608786505319974948498`35.}, {4, 0.1472610275223323382651759846329784927056430492700899336219`35.}, {7, 0.0522096126394425655340239042320925681430535814764949657287`35.}, {10, 0.0197215428156869480356026620646722197670633092026702825467`35.}, {16, 0.0030618217265287229459987187029086772676833771159536917165`35.}}}

## alpha_integral

<result query='integral exp(-7/(2 cos(x)^2))/arccos(sqrt(0.9)) from x=0 to arccos(sqrt(0.9))' url='https://www.wolframalpha.com/input?i=integral+exp%28-7%2F%282+cos%28x%29%5E2%29%29%2Farccos%28sqrt%280.9%29%29+from+x%3D0+to+arccos%28sqrt%280.9%29%29'>

# Definite integral

integral\_0^0.321751 exp(-7/(2 cos^2(x)))/(cos^(-1)(sqrt(0.9))) dx = 0.0267937

# Visual representation of the integral

![Visual representation of the integral](https://public6.wolframalpha.com/files/GIF_xwf6cr77le.gif)

</result>
