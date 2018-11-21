#!/usr/bin/python3

import turtle
from random import random
import math
import time

N = 5
turtle.speed(2)

data = [(random(),random()) for _ in range(N)]

def dist(p1,p2):
	return math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

EPSILON=0.0001
def getOrientation(p1, p2, p3):
	 val = (p2[1] - p1[1]) * (p3[0] - p2[0]) -\
		   (p2[0] - p1[0]) * (p3[1] - p2[0])
	 if (abs(val) < EPSILON):
		 return "COL"
	 elif (val > 0):
		 return "CW"
	 else:
		 return "CCW"

def circleRadius(p1, p2, p3):
	x1 = p1[0];  x2 = p2[0];  x3 = p3[0];
	y1 = p1[1];  y2 = p2[1];  y3 = p3[1];

	k1 = x1*x1 + y1*y1;
	k2 = x2*x2 + y2*y2;
	k3 = x3*x3 + y3*y3;
	k4 = y1-y2;
	k5 = y1-y3;
	k6 = y2-y3;
	k7 = x1-x2;
	k8 = x1-x3;
	k9 = x2-x3;
	k10 = x1*x1;
	k11 = x2*x2;
	k12 = x3*x3;

	mnom = (k12*(-k4)+k11*k5-(k10+k4*k5)*k6);
	mdenom = (-k4)*x3+x2*k5+x1*(-k6);
	m = 0.5*mnom/mdenom;

	nnom = k1*(-k9)+k2*k8+k3*(-k7);
	ndenom = y1*(-k9)+y2*k8+y3*(-k7);

	n = 0.5*nnom/ndenom;


	radius = dist((m,n), p2);
	return ((m,n),radius);

def checkInside(middle,radius,points):
	for p in points:
		if dist(middle,p) < radius-EPSILON:
			return True
	return False
	
def minimalBoundingCircle(e, points):
	minPoint = None
	minradius = float("inf")
	
	for p in points:
		if p == e[0] or p == e[1]:
			continue

		if (getOrientation(e[0], e[1], p) != "CCW"):
			continue;

		middle,radius = circleRadius(e[0], e[1], p);
		if (radius < minradius ):
			if checkInside(middle,radius,points):
				pass
#				continue
			minPoint = p;
			minradius = radius;
	return minPoint

def addToAel(ael,edge):
	print("AEL: {}, edge: {}".format(strAEL(ael),strEdge(edge)))
	if edge in ael:
		print("Terrible things are happening...")
	swapped = (edge[1],edge[0])
	print("Swapped: {}".format(strEdge(swapped)))
	if swapped in ael:
		print("Found, removing")
		swidx = ael.index(swapped)
		del ael[swidx]
		return
	print("Not found, adding")
	ael.append(edge)

def delaunay(points):
	dt = list()

	p1 = points[0];
	p2 = None
	
	mindist = float("inf")
	for p in points:
		if p == p1:
			continue
		if (dist(p1,p) < mindist):
			p2 = p;
			mindist = dist(p1,p);
		
	e = (p1,p2)
	p = minimalBoundingCircle(e, points);

	if not p:
		e = (e[1],e[0])
		p = minimalBoundingCircle(e, points);
		e2 = (p1,p);
		e3 = (p,p2);
	else:
		e2 = (p2,p);
		e3 = (p,p1);
	
	if not p:
		print("p not found")


	dt.append(e);
	dt.append(e2);
	dt.append(e3);

	turtle.pencolor(0,1,0)
	drawEdge(e)
	drawEdge(e2)
	drawEdge(e3)
	turtle.pencolor(0,0,0)
	
	ael = list()

	ael.append(e);
	ael.append(e2);
	ael.append(e3);

	while len(ael) > 0:
		print("Fronta: {}".format(len(ael)))
		e = ael.pop()
		print("Edge: {}".format(e))
		e = (e[1],e[0])
		turtle.pencolor(1,0,0)
		drawEdge(e)
		turtle.pencolor(0,0,0)
		time.sleep(5)

		p = minimalBoundingCircle(e, points);
			#System.out.format("p = %h\n",p);
			#System.out.println(e);
		if not p:
			continue
		e2 = (p,e[0])
		e3 = (e[1],p)

		dt.append(e);
		dt.append(e2);
		dt.append(e3);

		drawEdge(e2);
		drawEdge(e3);
		drawEdge(e);
		time.sleep(10);

		
		addToAel(ael, e3);
		addToAel(ael, e2);

	return dt;

def drawPoints(points):
	for (i,p) in enumerate(points):
		turtle.goto(p[0]*300,p[1]*300)
		turtle.dot()
		turtle.write(i)


def drawTriangulation(edges):
	for e in edges:
		drawEdge(e)

def strEdge(e):
	return "{}--{}".format(data.index(e[0]),data.index(e[1]))
def strAEL(ael):
	return ", ".join(map(strEdge,ael))

def drawEdge(e):
	turtle.st()
	((sx,sy),(ex,ey)) = e
	turtle.goto(sx*300,sy*300)
	turtle.pendown()
	turtle.goto(ex*300,ey*300)
	turtle.dot()
	turtle.penup()
	turtle.ht()

turtle.penup()
drawPoints(data)
edges = delaunay(data)
drawTriangulation(edges)
turtle.exitonclick()
