import numpy as np
import math

R_MAX=300#distributionFunctionの最大値
DELTA_R=1#ヒストグラムのビンの幅
X=1794#系のサイズ
Y=762#系のサイズ

def distributionFunction(xs,ys,boundary=False):
  if boundary:
    distribution_functions=[distributionFunctionParticleBoundary(xs,ys,x,y) for x,y in zip(xs,ys)]
    return np.mean( distribution_functions,axis=0 )
  else:
    distribution_functions=[distributionFunctionParticle(xs,ys,x,y) for x,y in zip(xs,ys)]
    return np.mean( distribution_functions,axis=0 )
  
def distributionFunctionParticle(xs,ys,x,y):
  rho=float(len(xs)/X/Y)
  drs=np.sqrt( (xs-x)*(xs-x)+(ys-y)*(ys-y) )
  ns,rs=np.histogram( drs,bins=int(R_MAX/DELTA_R),range=(0,R_MAX) )
  return ns[1:]/(2*math.pi*rs[1:-1]*rho*DELTA_R)#0で割らないためにスライス#binの方が要素が一個多い

def distributionFunctionParticleBoundary(xs,ys,x,y):
     rho=float(len(xs)/X/Y)
     
     drs=np.sqrt( (xs-x)*(xs-x)+(ys-y)*(ys-y) )
     ns,rs=np.histogram( drs,bins=int(R_MAX/DELTA_R),range=(0,R_MAX) )
     arcs=np.array( [arc(x,y,r) for r in rs[1:-1]] )#binの方が要素が一個多い
     return ns[1:]/(arcs*rho*DELTA_R)#0で割らないためにスライス
def getKakudosa(theta1,theta2):#theta2-theta1を計算
    theta1=theta1%(2*math.pi)
    theta2=theta2%(2*math.pi)
    delta_theta=theta2-theta1
    return np.where(np.abs(delta_theta)<math.pi, delta_theta, delta_theta-np.sign(delta_theta)*2*math.pi)*1.0
 def arc(x,y,r):
     def solve(cx, cy, r, line):#https://tjkendev.github.io/procon-library/python/geometry/circle_line_cross_point.html
         [(x1, y1), (x2, y2)]=line
         xd = x2 - x1; yd = y2 - y1
         X = x1 - cx; Y = y1 - cy
         a = xd**2 + yd**2
         b = xd * X + yd * Y
         c = X**2 + Y**2 - r**2
         # D = 0の時は1本で、D < 0の時は存在しない
         D = b**2 - a*c
         if D>0:
             s1 = (-b + math.sqrt(D)) / a
             s2 = (-b - math.sqrt(D)) / a
             return [[x1 + xd*s1, y1 + yd*s1], [x1 + xd*s2, y1 + yd*s2]]
         else:
             return []
     
     lines = [[(0, 0),(0, Y)],[(0, 0),(X,0)],[(X,0),(X,Y)],[(0,Y),(X,Y)]]
     pointss=[solve(x,y,r,line) for line in lines]
     
     flag=all( [all( [( 0<=point[0]<=X and 0<=point[1]<=Y) for point in points] ) for points in pointss] )
     if flag:
         theta=2*math.pi
         for points in pointss:
             if len(points)==2:
                 dr1=np.array(points[0])-np.array([x,y])
                 theta1=math.atan2(dr1[1],dr1[0])
                 dr2=np.array(points[1])-np.array([x,y])
                 theta2=math.atan2(dr2[1],dr2[0])
                 theta-=abs(getKakudosa(theta1,theta2))
         return theta*r
     else:
         points=pointss[0]+pointss[1]+pointss[2]+pointss[3]
         points=[point for point in points if 0<=point[0]<=X and 0<=point[1]<=Y]
         if len(points)==2:
             dr1=points[0]-np.array([x,y])
             theta1=math.atan2(dr1[1],dr1[0])
             dr2=points[1]-np.array([x,y])
             theta2=math.atan2(dr2[1],dr2[0])
             theta=abs(getKakudosa(theta1,theta2))
             return theta*r
         else:
             print("error...?",x,y,r,points)
             return 2*math.pi*r
