from logging import exception
from msilib.schema import Class
from operator import index
import os
import sys
from certifi import where
import ifcopenshell
import ifcopenshell.util
import ifcopenshell.mvd
import ifcopenshell.util.element   
import ifcopenshell.util.selector
import ifcopenshell.geom
import matplotlib
import matplotlib.pyplot as plt
import MySQLdb
import scipy
import alphashape
import scipy.spatial
import pandas as pd
import numpy as np
import math
import pylab
from shapely import *
import shapely.geometry
from descartes import PolygonPatch
from itertools import zip_longest as zl
from scipy.spatial import ConvexHull, convex_hull_plot_2d


def display(arg):
    for i in arg:
        print(i)
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ 
    https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def rotate_pt(xy, radians, origin=(0, 0)):
    """
    https://gist.github.com/LyleScott/e36e08bfb23b1f87af68c9051f985302
    Rotate a point around a given point.
    
    I call this the "high performance" version since we're caching some
    values that are needed >1 time. It's less readable than the previous
    function but it's faster.
    """
    x, y = xy
    offset_x, offset_y = origin
    adjusted_x = (x - offset_x)
    adjusted_y = (y - offset_y)
    cos_rad = math.cos(radians)
    sin_rad = math.sin(radians)
    qx = offset_x + cos_rad * adjusted_x + sin_rad * adjusted_y
    qy = offset_y + -sin_rad * adjusted_x + cos_rad * adjusted_y

    return qx, qy
def grouped(iterable, n):
    'https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list'
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zl(*[iter(iterable)]*n)

#def isBetween(a, b, c):
    crossproduct = (c.y - a.y) * (b.x - a.x) - (c.x - a.x) * (b.y - a.y)

    # compare versus epsilon for floating point values, or != 0 if using integers
    if abs(crossproduct) > epsilon:
        return False

    dotproduct = (c.x - a.x) * (b.x - a.x) + (c.y - a.y)*(b.y - a.y)
    if dotproduct < 0:
        return False

    squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y)
    if dotproduct > squaredlengthba:
        return False

    return True

class Wall:
    def __init__(CurWall,Guid,AxPoints,BoundPoints,tag):
        CurWall.id = Guid
        CurWall.AxPoints= AxPoints
        CurWall.BPoints=BoundPoints
        CurWall.Tag = tag

    def setExternalPts(CurWall,ExtPts):
        CurWall.Ext = ExtPts
    def setInternalPts(CurWall,IntPts):
        CurWall.Int = IntPts
NormX=(1,0,0)
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\SiteTest23.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\IFC_tt9.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\COL_4BEDROOM DUPLEX.ifc')
IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\LV_EdificioHabitacional_R07.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\ClosePointTest.ifc')
#display(IfcFile)

Walls= IfcFile.by_type('ifcWall')
print('**************************************walls****************************************************************')
#display (Walls)       
TestWall=Walls[0]
psets=ifcopenshell.util.element.get_psets(TestWall)
#print(psets)
ExternalWalls=[]
InternalWalls=[]
for i in Walls:
    PsetCom=ifcopenshell.util.element.get_psets(i)['Pset_WallCommon']
    if PsetCom['IsExternal'] == True:
        ExternalWalls.append(i)
    else:
        InternalWalls.append(i)
#display(ExternalWalls)
#print('***************************************Internal**********************************************************************************************')
#display(InternalWalls)
#print('Walls per Level************************************************************************************************************************')
WallSet=ExternalWalls
WallLev={}
for i in WallSet:
    if i.ContainedInStructure[0].RelatingStructure.LongName not in WallLev:
        WallLev[i.ContainedInStructure[0].RelatingStructure.LongName]=list()
        WallLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
    elif i.ContainedInStructure[0].RelatingStructure.LongName in WallLev:
        WallLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
    else:
        continue
#WallSet=ExternalWalls
WallSet=WallLev['1º PAVIMENTO N.A.']
#WallSet=WallLev['Level 1']
#print(WallLev)
#display(WallLev['Level 1'])
#print('xxxxxxxxxxxxxxxxxx************************xxxxxxxxxxxxxxxxxxx')
# display(WallLev['Level 2'])

AllPoints=[]
Vectors=[]
Angles=[]
WallAxPt=[]
ConWalls=[]
Location=[]
PtLst=[]
ItCount=0
WallEval=[]
for i in WallSet:
    WallEval.append(i)


# Check all the walls connected in certain DataSet, the conditions can be modified
WallLoop=[]
while len(WallEval)>0:
    TestWall=WallEval[0]
    x=0
    ConWalls=[]
    while TestWall not in ConWalls:
        if x==0:
            ConWalls.append(TestWall)
            WallEval.remove(TestWall)
            NextWall=[]
            print(TestWall in ConWalls)
            print('TestWall:', TestWall)
            for i in TestWall.ConnectedTo:
                print(x)
                print('this is i',i,'xxxxxxxxxxxxxxxxxxxxxxxxx')
                y=i.RelatedElement
                x+=1 
                if y not in ConWalls and y in WallSet:
                    print('Got into loop and this is y:',y)
                    TestWall=y
                    WallEval.remove(TestWall)
                    NextWall.append(y)
                    print('******New TestWall:', TestWall)
                    break
                else:                
                    continue
            if bool(NextWall) == False:
                try:
                    for i in TestWall.ConnectedFrom:
                        print(x)
                    print('this is i',i,'xxxxxxxxxxxxxxxxxxxxxxxxx')
                    y=i.RelatedElement
                    x+=1 
                    if y not in ConWalls and y in WallSet:
                        print('Got into loop and this is y:',y)
                        TestWall=y
                        WallEval.remove(TestWall)
                        NextWall.append(y)
                        print('******New TestWall:', TestWall)
                        break
                    else:                
                        continue
                except:
                    print("Wall doesn't seem to be connected")
        if x>0:
            print('********* loop of x>0***************')
            for i in NextWall:
                ConWalls.append(i)
            NextWall=[]
            print(TestWall in ConWalls)
            print('TestWall:', TestWall)
            for i in TestWall.ConnectedTo:
                print(x)
                print('this is i',i,'xxxxxxxxxxxxxxxxxxxxxxxxx')
                y=i.RelatedElement
                x=x+1 
                if y not in ConWalls and y in WallSet:
                    print('Got into loop and this is y:',y)
                    TestWall=y
                    WallEval.remove(TestWall)
                    NextWall.append(y)
                    print('******New TestWall:', TestWall)
                    break
                else:                
                    continue
            if bool(NextWall) == False:
                try:
                    for i in TestWall.ConnectedFrom:
                        print(x)
                        print('this is i',i,'xxxxxxxxxxxxxxxxxxxxxxxxx')
                        y=i.RelatingElement
                        x=x+1 
                        if y not in ConWalls and y in WallSet:
                            print('Got into loop and this is y:',y)
                            TestWall=y
                            WallEval.remove(TestWall)
                            NextWall.append(y)
                            print('******New TestWall:', TestWall)
                            break
                        else:                
                            continue
                    if bool(NextWall) == False:
                        print('No wall connected, changing TestWall')
                        TestWall=ConWalls[((len(ConWalls))-1)]
                except:
                    print("Wall doesn't seem to be connected")
                    continue 
    print('this are the connected walls+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    display(ConWalls)
    WallLoop.append(ConWalls)
    #Get the points of the profile of the wall checking if it is a sweepSolid (arbitratyProfile or RectProfile) and then checking if its a composite profile
WallObjects=[]
LoopWalls=[]
ExtShapes=[]
IntShapes=[]
for Loop in WallLoop:
    IntOrd=[]
    ExtOrd=[]
    AxPtList=[]
    WallAxPt=[]
    AllPoints=[]
    PtLst=[]
    if bool(WallObjects) != False:
        LoopWalls.append(WallObjects) 
    WallObjects=[]
    for i in Loop:
         
        PPoint=i.ObjectPlacement.RelativePlacement.Location.Coordinates
        Location.append(PPoint)
        try:
            ItCount=ItCount+1
            print(ItCount)
            PPoint=i.ObjectPlacement.RelativePlacement.Location.Coordinates
            print('--------------------------------------Este es el PPoint',PPoint)
            print(i)
            print(i.ObjectPlacement)
            print(i.ObjectPlacement.RelativePlacement)
            try:
                directionX=i.ObjectPlacement.RelativePlacement.RefDirection.DirectionRatios
            except:
                print('No axis or refdirection found, using default placement')
                directionX=(1,0,0)
            Angle=angle_between(NormX,directionX)
            Angles.append(Angle)
            AxisPts=i.Representation.Representations[0].Items[0].Points
            print('esta es la dirección {}, el ángulo {} y los puntos {}'.format(directionX,Angle,AxisPts))
            print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',AxisPts)
            AxTempPt=[]
            for j in AxisPts:
                AxCoord=j.Coordinates
                
                if directionX[1] >= 0:
                    RAxCoord=rotate_pt(AxCoord,-Angle)
                    #print('22222222222222222222222222222222222-Rotated  Axcoords',RAxCoord)
                    AxConvPt=(RAxCoord[0]+PPoint[0],RAxCoord[1]+PPoint[1])
                    ##print('Coordenadas de RelPlace',PPoint)
                    AxTempPt.append(AxConvPt)
                elif directionX[1] <0:
                    RAxCoord=rotate_pt(AxCoord,Angle)
                    #print('22222222222222222222222222222222222-Rotated Ax coords',RAxCoord)
                    AxConvPt=(RAxCoord[0]+PPoint[0],RAxCoord[1]+PPoint[1])
                    ##print('Coordenadas de RelPlace',PPoint)
                    AxTempPt.append(AxConvPt)
                #print('3333333333333333333333333333333333333-Normalized Ax Coords', AxConvPt)
                WallAxPt.append(AxTempPt)
        except Exception as e:
            print(e)
            print('hubo un error con el muro {}'.format(i.Tag))
            pass
        try:
            try:
                directionX=i.ObjectPlacement.RelativePlacement.RefDirection.DirectionRatios
            except:
                print('No axis or refdirection found, using default placement')
                directionX=(1,0,0)
            Vectors.append(directionX)
            Angle=angle_between(NormX,directionX)
            Angles.append(Angle)
            BTempPt=[]
            print('*******************************************************************this is the vector',directionX, 'and the angle',Angle,' or ',math.degrees(Angle))
            if i.Representation.Representations[1].RepresentationType == 'SweptSolid':
                print('SweptSolid with rectangular shape found')
                try:
                    DimX=i.Representation.Representations[1].Items[0].SweptArea.XDim
                    DimY=i.Representation.Representations[1].Items[0].SweptArea.YDim
                    Points=[(0,-(DimY/2)),(DimX,-(DimY/2)), (0,DimY/2),(DimX,DimY/2)]
                    for j in Points:    
                        Coord=j
                        #print('11111111111111111111111111111111111111111111111111 these are the relative coords',Coord)
                        if directionX[1] >= 0:
                            RCoord=rotate_pt(Coord,-Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            ##print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        elif directionX[1] <0:
                            RCoord=rotate_pt(Coord,Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            ##print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        #print('3333333333333333333333333333333333333-Normalized Coords', ConvPt)
                
                        print('**********************************xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx****************************************')
                    AllPoints.append(BTempPt)
                    print('these are the points of',i,'**********************************************************************')
                    display(BTempPt)
                    print('end of the points of {}'.format(i))
                except:
                    print('No Rectangular Profile')
                    Points=i.Representation.Representations[1].Items[0].SweptArea.OuterCurve.Points
                    for j in Points:
                        Coord=j.Coordinates
                        #print('11111111111111111111111111111111111111111111111111 these are the relative coords',Coord)
                        if directionX[1] >= 0:
                            RCoord=rotate_pt(Coord,-Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            #print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        elif directionX[1] <0:
                            RCoord=rotate_pt(Coord,Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            #print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        #print('3333333333333333333333333333333333333-Normalized Coords', ConvPt)
                
                        print('**********************************xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx****************************************')
                    AllPoints.append(BTempPt)
                    print('these are the points of',i,'**********************************************************************')
                    display(BTempPt)
                    print('end of the points of {}*******************************************************************'.format(i))
            else:
                try:
                    Points=i.Representation.Representations[1].Items[0].FirstOperand.SweptArea.OuterCurve.Points
                    for j in Points:
                        Coord=j.Coordinates
                        #print('11111111111111111111111111111111111111111111111111 these are the relative coords',Coord)
                        if directionX[1] >= 0:
                            RCoord=rotate_pt(Coord,-Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            ##print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        elif directionX[1] <0:
                            RCoord=rotate_pt(Coord,Angle)
                            #print('22222222222222222222222222222222222-Rotated coords',RCoord)
                            ConvPt=(RCoord[0]+PPoint[0],RCoord[1]+PPoint[1])
                            ##print('Coordenadas de RelPlace',PPoint)
                            BTempPt.append(ConvPt)
                        #print('3333333333333333333333333333333333333-Normalized Coords', ConvPt)
                
                        print('**********************************xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx****************************************')
                    AllPoints.append(BTempPt)
                    print('these are the points of',i,'**********************************************************************')
                    display(BTempPt)
                    print('end of the points of {}'.format(i))
                except:
                    print('No Clipping Boundary found')
        
        except:
            pass
        BoundPoint=[]
        for x in BTempPt:
            tup=(round(x[0],6)),round(x[1],6)
            BoundPoint.append(tup)
        CWall=Wall(i.GlobalId,AxTempPt,BoundPoint,i.Tag)
        WallObjects.append(CWall)
    print('these are the wall objects****************************************************')
    for i in WallObjects:
        print('this is the wall {}'.format(i.Tag))
        print (i.BPoints)
      
    print('these are the raw points!!!!!!!!**********************************************************************************')
    display(AllPoints)
    print('end of raw points**********************************************************************************************')
    for i in AllPoints:
        for j in i:
            PtLst.append(j)
    AxPtList=[]
    for i in WallAxPt:
        for j in i:
            if j not in AxPtList:
                AxPtList.append(j) 

    print('xxxxxxxxxxxxxxxxxxxxxx**********************************************xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')

    display(PtLst)
    print('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')
    display(Vectors)
    display(Angles)
    Degrees=[]
    for i in Angles:
        Deg=math.degrees(i)
        Degrees.append(Deg) 
        #print('Loc Element {}Ángulo rad={} Ángulo Grad={} Vector={}'.format(Location[(Angles.index(i))],i,Deg,Vectors[(Angles.index(i))]))
    display(Location)    

    display(AxPtList)
    PolyPoints=[]
    for i in Loop:
        PolyPoints.append(i.ObjectPlacement.RelativePlacement.Location.Coordinates)
    print('These are the polypoints**********************')
    print(PtLst)
    OrdPts=[]
    #This iteration evaluate two connected walls by its axis points, then determines the order inside the wall to perform the connection.
    for i,j,l,m in zip(AxPtList[0::4],AxPtList[1::4],AxPtList[2::4],AxPtList[3::4]):
        print('these are the points**************************************',i,j,l,m)
        dist1=math.dist(i,l)
        dist2=math.dist(i,m)
        dist3=math.dist(j,l)
        dist4=math.dist(j,m)
        if dist1 < dist2 and dist1 <dist3 and dist1< dist4:
            OrdPts.append(j)
            OrdPts.append(i)
            OrdPts.append(l)
            OrdPts.append(m)
        elif dist2 < dist1 and dist2 <dist3 and dist2< dist4:
            OrdPts.append(j)
            OrdPts.append(i)
            OrdPts.append(m)
            OrdPts.append(l)
        elif dist3 < dist1 and dist3 <dist2 and dist3< dist4:
            OrdPts.append(i)
            OrdPts.append(j)
            OrdPts.append(l)
            OrdPts.append(m)
        elif dist4 < dist1 and dist4 <dist2 and dist4< dist3:
            OrdPts.append(i)
            OrdPts.append(j)
            OrdPts.append(m)
            OrdPts.append(l)
    for i, j, l, m in grouped(AxPtList,4):
        temp=[]
        print(i,j,l,m)
        if m is None:        
            try:
                dist1=math.dist(i,l)
                dist2=math.dist(j,l)
                if dist1 < dist2:
                    OrdPts.append(j)
                    OrdPts.append(i)
                    OrdPts.append(l)  
                else:
                    OrdPts.append(i)
                    OrdPts.append(j)
                    OrdPts.append(l)       
            except:
                pass
        if l is None:
            try:
                tpPt=OrdPts[0]
                dist1=math.dist(i,tpPt)
                dist2=math.dist(j,tpPt)
                if dist1<dist2:
                    OrdPts.append(j)
                    OrdPts.append(i)
                    break
                else:
                    OrdPts.append(i)
                    OrdPts.append(j)
                    break
            except:
                OrdPts.append(i)
                OrdPts.append(j)
                pass
        elif j is None:
            OrdPts.append(i)
                    
            
            
        
                
    print('These are the ordinated pointsxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')






    fig, ax = plt.subplots()



    #ax.scatter(*zip(*PtLst),c='blue')
    ax.scatter(*zip(*OrdPts),c='red')
        



    #Find the centroid of the polygon: https://stackoverflow.com/questions/10846431/ordering-shuffled-points-that-can-be-joined-to-form-a-polygon-in-python    
    cent=(sum([p[0] for p in OrdPts])/len(OrdPts),sum([p[1] for p in OrdPts])/len(OrdPts))
    print(cent)
    #AxPtList.sort(key=lambda p: math.atan2(p[1]-cent[1],p[0]-cent[0]))

    Axs= shapely.geometry.Polygon(OrdPts)
    Orig=shapely.geometry.Polygon(PtLst)
    ExtPts=[]   
    IntPts=[]

    for i in WallObjects:
        WallExt=[]
        WallInt=[]
        for j in i.BPoints:
            tempPoint=shapely.geometry.Point(j)
            if Axs.contains(tempPoint):
                IntPts.append(j)
                WallInt.append(j)
            else:
                ExtPts.append(j)
                WallExt.append(j)
        i.setExternalPts(WallExt)
        i.setInternalPts(WallInt)
        print(i.Ext)        
                

    FirLst=[]
    TmpLst=[]
    for i in ExtPts:
        tup=(round(i[0],6)),round(i[1],6)
        if tup not in TmpLst:
            TmpLst.append(tup)

       
    print('these are the temp points******************************************')
    display(TmpLst)
    ExtOrd=[]
    grpset=0
    lastPt=()


    # Connect the exterior points based on a tuple of points
    addpt=0
    numwall=0
    group=0
    iteration=0

    for x,y in zip(WallObjects[0::2],WallObjects[1::2]):
        if iteration==0:
            dist1=math.dist(x.Ext[0],y.Ext[0])
            dist2=math.dist(x.Ext[1],y.Ext[0])
            dist3=math.dist(x.Ext[0],y.Ext[-1])
            dist4=math.dist(x.Ext[1],y.Ext[-1])
            if dist1 < dist2 and dist3 < dist4:
                lastPt=x.Ext[1]
            else:
                lastPt=x.Ext[0]
            iteration=iteration+1
        if iteration >0:
            for i,j in zip(TmpLst[0::2],TmpLst[1::2]):
                if bool(lastPt) == True:
                    break
                failcount=0
                SetPt=[i,j]
                if i in x.BPoints:
                    mur=1
                elif i in y.BPoints:
                    mur=2
                else:
                    print('x Point is not in the pair of walls')
                    failcount=1
                    pass
                if j in y.BPoints:
                    mur2=2
                elif j in x.BPoints:
                    mur2=1
                else:
                    failcount=failcount+1
                if failcount == 2:
                    break
                print('these are the points', SetPt,'in group: ',grpset,'....el punto {} pertenece al muro {} y el punto {} pertenece al muro {}'.format(i,mur,j,mur2))
                grpset=grpset+1
                dists=[]
                distsW1=[]
                distsW2=[]
                if lastPt == ():            
                    dist1=math.dist(i,TmpLst[2])
                    dist2=math.dist(i,TmpLst[3])
                    dist3=math.dist(j,TmpLst[2])
                    dist4=math.dist(j,TmpLst[3])
                    if dist1<=dist3 :
                        lastPt=j
                        
                        if dist4< dist2:
                            lastPt=i
                            
                    elif dist3<dist1 :
                        lastPt=i
                
        evalpts=[]
        for i in x.Ext:
            if i in TmpLst:
                evalpts.append(i)
        for i in y.Ext:
            if i in TmpLst:
                evalpts.append(i)
        print('these are the evalpoints:{}'.format(evalpts))
        dists=[]
        while len(evalpts) >0:
            print('while iteration started')
    
            indexes=[]
            for i in enumerate(evalpts):
                if i[1] in x.BPoints:
                    dist=math.dist(i[1],lastPt) 
                    dists.append(dist)
                    indexes.append(i[0])  
            if bool(dists) == False:
                for i in enumerate(evalpts):
                    if i[1] in y.BPoints:
                        dist=math.dist(i[1],lastPt) 
                        dists.append(dist)
                        indexes.append(i[0])
                        ###########Revisar
            min_dist=min(dists)                
            indMin=dists.index(min_dist)
            ptIndx=indexes[indMin]
            ClPt=evalpts[ptIndx]
            lastPt=ClPt
            if lastPt not in ExtOrd:         
                ExtOrd.append(lastPt)
                print('aqui se añadió el punto :{} con coordenadas {}'.format(addpt,lastPt))
                addpt=addpt+1
            evalpts.remove(ClPt)
            dists=[]
    if len(WallObjects)%2>0:
        evalpts=[]
        for i in WallObjects[-1].BPoints:
            if i in TmpLst:
                evalpts.append(i)
        dists=[]
        while len(evalpts) >0:
            indexes=[]
            for i in enumerate(evalpts):
                dist=math.dist(i[1],lastPt) 
                dists.append(dist)
                indexes.append(i[0])  
            min_dist=min(dists)                
            indMin=dists.index(min_dist)
            ptIndx=indexes[indMin]
            ClPt=evalpts[ptIndx]
            lastPt=ClPt
            if lastPt not in ExtOrd:         
                ExtOrd.append(lastPt)
                print('aqui se añadió el punto :{} con coordenadas {}'.format(addpt,lastPt))
                addpt=addpt+1
            evalpts.remove(ClPt)
            dists=[]
    LastOrdPT=ExtOrd[0]
    '''for i in ExtOrd:
        if ExtOrd.index(i)==0:
            continue
        else:
            dist=math.dist(LastOrdPT,i)
            if dist <0.5:
                Ind=ExtOrd.index(LastOrdPT)
                ExtOrd.remove(LastOrdPT)
                LastOrdPT=i
            else:
                LastOrdPT=i'''  


    #ax.scatter(*zip(Axs.exterior.xy),c='red')
    print('These are the exterior points***********************************************************************')
    display(ExtPts)
    print('These are the exterior ordered points***********************************************************************')
    display(ExtOrd)
    ExtShp=shapely.geometry.Polygon(ExtOrd)
    ExtShapes.append(ExtShp)
    #Connect the interior points
    numwall=0
    TmpLst=[]
    for i in IntPts:
        tup=(round(i[0],6)),round(i[1],6)
        if tup not in TmpLst:
            TmpLst.append(tup)
    IntOrd=[]
    grpset=0
    lastPt=()
    iteration=0
    addpt=0
    for x,y in zip(WallObjects[0::2],WallObjects[1::2]):
        for i,j in zip(TmpLst[0::2],TmpLst[1::2]):
            if bool(lastPt) == True:
                break
            failcount=0
            SetPt=[i,j]
            if i in x.BPoints:
                mur=1
            else:
                mur=2
            if j in y.BPoints:
                mur2=2
            else:
                mur2=1
            print('these are the points', SetPt,'in group: ',grpset,'....el punto {} pertenece al muro {} y el punto {} pertenece al muro {}'.format(i,mur,j,mur2))
            grpset=grpset+1
            dists=[]
            distsW1=[]
            distsW2=[]
            
            if lastPt == ():
                dist1=math.dist(i,ExtOrd[0])
                dist2=math.dist(j,ExtOrd[0])
                if dist1<dist2:
                    lastPt=i
                else:
                    lastPt=j
        evalpts=[]
        for i in x.BPoints:
            if i in TmpLst:
                evalpts.append(i)
        for i in y.BPoints:
            if i in TmpLst:
                evalpts.append(i)
        print('these are the evalpoints:{}'.format(evalpts))
        dists=[]
        while len(evalpts) >0:

            print('while iteration started')

            indexes=[]
            for i in enumerate(evalpts):
                if i[1] in x.BPoints:
                    dist=math.dist(i[1],lastPt) 
                    dists.append(dist)
                    indexes.append(i[0])  
            if bool(dists) == False:
                for i in enumerate(evalpts):
                    if i[1] in y.BPoints:
                        dist=math.dist(i[1],lastPt) 
                        dists.append(dist)
                        indexes.append(i[0])
            min_dist=min(dists)                
            indMin=dists.index(min_dist)
            ptIndx=indexes[indMin]
            ClPt=evalpts[ptIndx]
            lastPt=ClPt
            if lastPt not in IntOrd:         
                IntOrd.append(lastPt)
                print('aqui se añadió el punto :{} con coordenadas {}'.format(addpt,lastPt))
                addpt=addpt+1
            evalpts.remove(ClPt)
            dists=[]
    print('This are the Interior points******************************************************')
    display(IntPts)
    IntShp=shapely.geometry.Polygon(IntOrd)
    IntShapes.append(IntShp)
    fig, ax = plt.subplots()
    ax.scatter(*zip(*ExtOrd),c='green')
    
    ax.scatter(*zip(*TmpLst),c='cyan')

    #ax.scatter(*zip(*IntOrd),c='black')
    print(Axs.wkt)
    print(Axs.is_valid)

    for i in enumerate(ExtOrd):
        ax.annotate(i[0],(i[1][0],i[1][1]))
        #ax.annotate(i[1],(i[1][0],i[1][1]))
    print(Axs.area)
    print('this is the floor area',ExtShp.area)
    ax.add_patch(PolygonPatch(Axs, fc='red', alpha=0.2))

    ax.add_patch(PolygonPatch(ExtShp, fc='blue', alpha=0.2))

    ax.add_patch(PolygonPatch(IntShp, fc='green', alpha=0.2))
    plt.show()

fig, ax = plt.subplots()     
ax.scatter(cent[0],cent[1],c='cyan')
ax.add_patch(PolygonPatch(Axs, fc='red', alpha=0.2))
for i in ExtShapes:
    ax.add_patch(PolygonPatch(i, fc='blue', alpha=0.2))
for i in IntShapes:
    ax.add_patch(PolygonPatch(i, fc='green', alpha=0.2))

plt.show()
