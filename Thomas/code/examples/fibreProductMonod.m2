--compute monodromy of fibre product.
restart
loadPackage("Monodromy")

R = CC[a_0..a_9][x_1,x_2,y_1,y_2]
F = {a_0+a_1*x_2+a_2*x_1+a_3*x_1*x_2+a_4*x_1^2,a_5+a_6*x_2+a_7*x_1+a_8*x_1*x_2+a_9*x_1^2,a_0+a_1*y_2+a_2*y_1+a_3*y_1*y_2+a_4*y_1^2,a_5+a_6*y_2+a_7*y_1+a_8*y_1*y_2+a_9*y_1^2}
baseP = apply(10,i->5*(random(CC)-random(CC)))

S = CC[t_0..t_3]
phi = map(S,R,{t_0,t_1,t_2,t_3}|baseP)
baseS = solveSystem(F/phi)/coordinates

bigD = select(baseS,s->abs(s#0-s#2)+abs(s#1-s#3)<1e-5)
baseS = select(baseS,s->not member(s,bigD))

M = monodromy(F,baseP,baseS)
monodromyLoop(M,50)
printCycleTally M
