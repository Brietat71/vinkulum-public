"""Contre-épreuve rationnelle de l'identité d'attribution discrète.

Cette identité ne suffit pas à prouver la qualité d'une trajectoire : les
résidus locaux peuvent être arbitrairement grands. Ils devront être bornés.
"""
from fractions import Fraction as F


def vec(*xs):return list(map(F,xs))
def add(a,b):return [x+y for x,y in zip(a,b)]
def sub(a,b):return [x-y for x,y in zip(a,b)]
def scale(s,a):return [s*x for x in a]
def eye():return [[F(i==j) for j in range(3)] for i in range(3)]
def mv(a,x):return [sum(u*v for u,v in zip(row,x)) for row in a]
def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def ma(a,b):return [add(x,y) for x,y in zip(a,b)]
def ms(s,a):return [scale(s,row) for row in a]
def cross(x,y):return [x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]
def skew(x):return [vec(0,-x[2],x[1]),vec(x[2],0,-x[0]),vec(-x[1],x[0],0)]


def inv(a):
    out=[list(row)+e for row,e in zip(a,eye())]
    for j in range(3):
        pivot=next(i for i in range(j,3) if out[i][j])
        out[j],out[pivot]=out[pivot],out[j]
        q=out[j][j];out[j]=[x/q for x in out[j]]
        for i in range(3):
            if i!=j:
                q=out[i][j];out[i]=sub(out[i],scale(q,out[j]))
    return [row[3:] for row in out]


class Champ:
    def __init__(self,inertie,frame):
        self.ji=inv(inertie);self.c=frame;self.ci=inv(frame)

    def f(self,x):
        p=mv(self.ci,x)
        return mv(self.c,cross(p,mv(self.ji,p)))

    def df(self,x):
        p=mv(self.ci,x)
        return mm(mm(self.c,ma(ms(-1,skew(mv(self.ji,p))),mm(skew(p),self.ji))),self.ci)


def verifie(a,b,x0,x1,y0,y1,h):
    mx=scale(F(1,2),add(x0,x1));my=scale(F(1,2),add(y0,y1))
    middle=scale(F(1,2),add(mx,my))
    jac=a.df(middle)
    assert sub(a.f(my),a.f(mx))==mv(jac,sub(my,mx))
    defect_a=sub(sub(x1,x0),scale(h,a.f(mx)))
    defect_b=sub(sub(y1,y0),scale(h,b.f(my)))
    left=ma(eye(),ms(-h/2,jac));right=ma(eye(),ms(h/2,jac))
    initial=sub(y0,x0);final=sub(y1,x1)
    parameter=scale(h,sub(b.f(my),a.f(my)))
    defect=sub(defect_b,defect_a)
    assert mv(left,final)==add(mv(right,initial),add(parameter,defect))
    li=inv(left)
    pieces=[mv(mm(li,right),initial),mv(li,parameter),mv(li,defect)]
    assert add(add(pieces[0],pieces[1]),pieces[2])==final
    return pieces


if __name__=='__main__':
    j=[vec('0.001',0,0),vec(0,'0.002',0),vec(0,0,'0.003')]
    # Transformation volontairement NON orthogonale, pour ne pas supposer
    # que les matrices stockées en binary64 sont des rotations exactes.
    c=[vec(1,F(1,7),0),vec(0,1,F(1,11)),vec(F(1,13),0,1)]
    a=Champ(j,eye());b=Champ(j,c)
    for k in range(1,21):
        x0=vec(F(k,1000),F(1,25),F(1,11000));x1=add(x0,vec(F(1,10**8),F(k,10**9),F(-1,10**9)))
        y0=add(x0,vec(F(1,10**16),F(1,10**17),F(1,10**18)))
        y1=add(x1,vec(F(1,10**12),F(-1,10**14),F(1,10**13)))
        verifie(a,b,x0,x1,y0,y1,F(1,1000))
    print('20 identités exactes, transformations non orthogonales incluses ; pas de certificat de trajectoire native.')
