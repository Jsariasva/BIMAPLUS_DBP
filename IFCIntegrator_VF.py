from cgi import test
from logging import exception
from msilib.schema import Class
from operator import index
import os
from re import I
import sys
#from certifi import where
import ifcopenshell
import ifcopenshell.util
import ifcopenshell.mvd
import ifcopenshell.util.element   
import ifcopenshell.util.selector
import ifcopenshell.util.placement
import ifcopenshell.util.representation
import ifcopenshell.geom
#from ifcopenshell import ids as ids
import matplotlib
import matplotlib.pyplot as plt
import MySQLdb
import scipy
import alphashape
import scipy.spatial
import pandas as pd
import numpy as np
import math
import OCC
import pylab
from shapely import *
import shapely.geometry
from shapely.ops import cascaded_union
from shapely.ops import unary_union
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

def Area(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

class Slab:
    def __init__(CurrSlab,Guid,Area,BoundPoints,tag,Shp,VoidPts=None,Voids=None):
        CurrSlab.id = Guid
        CurrSlab.Area= Area
        CurrSlab.BPoints=BoundPoints
        CurrSlab.Tag = tag
        CurrSlab.Shp = Shp
        CurrSlab.VoidPts = VoidPts
        CurrSlab.Voids = Voids

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
AreaDoLote=0
AreaDeImpl=0
AreaDeLogra=0
AreaImp=0
AreaBrutaC=0
LevABC={}
ImpShp=[]

NormX=(1,0,0)
Location=[]
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\SiteTest17.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\IFC_tt7.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\COL_4BEDROOM DUPLEX.ifc')
IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\LV_EdificioHabitacional_R10.ifc')
SelLev=[]
Levels=IfcFile.by_type('IfcBuildingStorey')
ChsnLvl=[]
while SelLev != 'N':
    print('These are the levels of the project:')
    for i in Levels:
        print(i.LongName)
    print('Please select a Level from the list above or type "ALL" for all the levels listed')
    
    SelLev=input('')
    if SelLev == "ALL":
        for i in Levels:
            Levels.remove(i)
            ChsnLvl.append(i)
    else:
        for i in Levels:
            if SelLev == i.LongName:
                Levels.remove(i)
                ChsnLvl.append(i)

    
    print('These are the selected levels')
    for i in ChsnLvl:
        print(i.LongName)
    print('Do you want to select more? (Y/N)')
    SelLev=input('')

LvlSlabShp={}
LvlWallExtShp={}
LvlWallIntShp={}
LvlImpShp={}


#******************************************This is the Site Module*******************************************************
fig, ax = plt.subplots()
sites = IfcFile.by_type('IfcSite')
display(sites)
SiteCoord=sites[0].Representation.Representations[1].Items
print(SiteCoord)
SitePoints=[]
for i in SiteCoord:
    a=[]
    points=i.Points
    for j in points:
        coord=j.Coordinates
        a.append(coord)
    SitePoints.append(a)
#******************************Graphic Display********************************************************
print('**********************these are the points****************************************************************************')
display(SitePoints)
for i in SitePoints:
    count=0
    for j in i:
        ax.scatter(j[0],j[1],c='green')
        ax.annotate(count,(j[0],j[1]),fontsize=8)
        count=count+1
print('**********************This is the Area********************************************************************************')
PlotShps=[]
for i in SitePoints:
    PlotArea=Area(i)
    PlotShp=shapely.geometry.Polygon(i)
    ax.add_patch(PolygonPatch(PlotShp, fc='green', alpha=0.2))
    print(PlotArea)
    PlotShps.append(PlotShp)
AreaDoLote=PlotArea    
plt.show()

for Lev in ChsnLvl:
    #****************************************This is the Slab Module!****************************************************************
    #**************************************Filter all the slabs by shape, including curved profiles and BREPS************************
    fig, ax = plt.subplots()
    Slabs=IfcFile.by_type('IfcSlab')
    CurvedSlabs=[]
    OrtoSlabs=[]
    Other=[]
    for i in Slabs:
        try:
            print(i.Representation.Representations[0].Items[0].SweptArea)
        except:
            pass
        if i.Representation.Representations[0].RepresentationType=='Brep':
            print(i)
            print('Brep found************************************************************************************')

            continue
        elif i.Representation.Representations[0].Items[0].SweptArea.is_a('IfcRectangleProfileDef'):
            print(i)
            print('Rectangular Profile found') 
            OrtoSlabs.append(i)
            continue
        elif i.Representation.Representations[0].Items[0].SweptArea.is_a('IfcArbitraryClosedProfileDef'):
            try:
                BoundType=i.Representation.Representations[0].Items[0].SweptArea.OuterCurve
                print(BoundType.is_a('IfcCompositeCurve'))
                if BoundType.is_a('IfcCompositeCurve'):
                    CurvedSlabs.append(i)
                elif BoundType.is_a('IfcPolyline'):
                    OrtoSlabs.append(i)
                continue
            except:
                continue
        else:
            Other.append(i)
            continue
    print('these are the curved slabs',CurvedSlabs)
    print('these are the ortoSlabs', OrtoSlabs)
    #******************************************Define a dictionary including all the ORTO entities per host level****************************************
    SlabLev={}
    for i in OrtoSlabs:
        try:
            if i.ContainedInStructure[0].RelatingStructure.LongName not in SlabLev:
                SlabLev[i.ContainedInStructure[0].RelatingStructure.LongName]=list()
                SlabLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
            elif i.ContainedInStructure[0].RelatingStructure.LongName in SlabLev:
                SlabLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
            else:
                continue
        except Exception as e:
            print(e)
            if i.Decomposes[0].RelatingObject.ContainedInStructure[0].RelatingStructure.LongName not in SlabLev:
                SlabLev[i.Decomposes[0].RelatingObject.ContainedInStructure[0].RelatingStructure.LongName]=list()
                SlabLev[i.Decomposes[0].RelatingObject.ContainedInStructure[0].RelatingStructure.LongName].append(i)
            elif i.Decomposes[0].RelatingObject.ContainedInStructure[0].RelatingStructure.LongName in SlabLev:
                SlabLev[i.Decomposes[0].RelatingObject.ContainedInStructure[0].RelatingStructure.LongName].append(i)


            
    #testSlab=OrtoSlabs[0]
    if bool(Lev)!= False:
        SlabSet=SlabLev[Lev.LongName]
    else:
        SlabSet=OrtoSlabs

    print(SlabSet)
    #*************************Create the shapes based on the geometry of the IFC Entity and include them into a list**************************************
    SlabEn=[]
    Shapes=[]
    for i in SlabSet:
        print(i.Representation.Representations[0].Items[0].SweptArea)

    for i in SlabSet:
        SlabPoints=[]
        VoidPts=[]
        SlabLoc=i.ObjectPlacement.RelativePlacement.Location.Coordinates
        if i.PredefinedType=='ROOF':
            SlabLoc=i.Decomposes[0].RelatingObject.ObjectPlacement.RelativePlacement.Location.Coordinates
        if i.Representation.Representations[0].Items[0].SweptArea.is_a('IfcArbitraryClosedProfileDef'):
            print(i)
            SlbPts=i.Representation.Representations[0].Items[0].SweptArea.OuterCurve.Points
            print(SlbPts)
            SlbPts2=[]        
            PlacePoint=i.Representation.Representations[0].Items[0].Position.Location.Coordinates
            print(PlacePoint)
            Location.append(PlacePoint)
            
            for j in SlbPts:
                pts=j.Coordinates

                try:
                    print(i.Representation.Representations[0].Items[0].Position.Axis[0])
                except:
                    pass
                try:
                    if i.Representation.Representations[0].Items[0].Position.Axis[0][2] < 0:
                        print('se evaluó desde abajo')
                        pts=((j.Coordinates[0]),(-j.Coordinates[1]))
                        SlbPts2.append(pts)
                    else:
                        print('se evaluó desde arriba')
                        pts=((j.Coordinates[0]),(j.Coordinates[1]))
                        SlbPts2.append(pts)
                except:
                    pts=((j.Coordinates[0]),(j.Coordinates[1]))
                    SlbPts2.append(pts)
            print('*********************************These are the SlbPts2:',SlbPts2)
            try:
                directionX=i.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios
            except:
                print('No axis or refdirection found, using default placement')
                directionX=(1,0,0)
            try:
                if i.Representation.Representations[0].Items[0].Position.Axis[0][2] < 0:
                    Angle=-(angle_between(NormX,directionX))
                else:
                    Angle=(angle_between(NormX,directionX))
            except:
                Angle=-(angle_between(NormX,directionX))
                pass
            Deg=math.degrees(Angle)
            print ('this is the angle in degrees {}'.format(Deg))
            print('this is the angle {}'.format(Angle))
            display(SlbPts)

            for k in SlbPts2:
                Coord=k
                if directionX[1] >= 0:
                    print('pt before rotation {}'.format(k))
                    Ri=rotate_pt(Coord,Angle)
                    print('pt after rotation {}'.format(Ri))
                    RealPt=(Ri[0]+PlacePoint[0]+SlabLoc[0],Ri[1]+PlacePoint[1]+SlabLoc[1])
                    SlabPoints.append(RealPt)
                elif directionX[1] <0:
                    Ri=rotate_pt(Coord,Angle)
                    RealPt=(Ri[0]+PlacePoint[0]+SlabLoc[0],Ri[1]+PlacePoint[1]+SlabLoc[1])
                    SlabPoints.append(RealPt)
        elif i.Representation.Representations[0].Items[0].SweptArea.is_a('IfcRectangleProfileDef'):
            Shape=i.Representation.Representations[0].Items[0].SweptArea
            SlbPts=(-(Shape.XDim)/2,-(Shape.YDim)/2),(-(Shape.XDim)/2,(Shape.YDim)/2),((Shape.XDim)/2,(Shape.YDim)/2),((Shape.XDim)/2,-(Shape.YDim)/2)

            SlbPts2=[]        
            PlacePoint=i.Representation.Representations[0].Items[0].Position.Location.Coordinates

            try:
                directionX=i.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios
                if i.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1] != 0:
                    directionX=i.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios
                    directionX=(directionX[0],directionX[1],0)
            except:
                print('No axis or refdirection found, using default placement')
                directionX=(1,0,0)
            
            try:
                if i.Representation.Representations[0].Items[0].Position.Axis[0][2] < 0:
                    print('se evaluó desde abajo')
                    Angle=-(angle_between(NormX,directionX))
                    print('This is the rotation angle with negative axis***************',Angle)
                else:
                    Angle=(angle_between(NormX,directionX))
                    print('This is the rotation angle***************',Angle)
            except:
                Angle=(angle_between(NormX,directionX))
                print('This is the rotation angle in exception***************',Angle)
                pass
            for j in SlbPts:
                try:
                    print(i.Representation.Representations[0].Items[0].Position.Axis[0])
                except:
                    pass
                try:
                    if i.Representation.Representations[0].Items[0].Position.Axis[0][2] < 0:
                        print('se evaluó desde abajo')
                        pts=((j[0]),(-j[1]))
                        SlbPts2.append(pts)
                    else:
                        print('se evaluó desde arriba')
                        pts=((j[0]),(j[1]))
                        SlbPts2.append(pts)
                except:
                    pts=((j[0]),(j[1]))
                    SlbPts2.append(pts)
            for k in SlbPts2:
                Coord=k
                if directionX[1] >= 0:
                    print('pt before rotation {}'.format(k))
                    Ri=rotate_pt(Coord,Angle)
                    print('pt after rotation {}'.format(Ri))
                    RealPt=(Ri[0]+PlacePoint[0]+SlabLoc[0],Ri[1]+PlacePoint[1]+SlabLoc[1])
                    SlabPoints.append(RealPt)
                elif directionX[1] <0:
                    Ri=rotate_pt(Coord,Angle)
                    RealPt=(Ri[0]+PlacePoint[0]+SlabLoc[0],Ri[1]+PlacePoint[1]+SlabLoc[1])
                    SlabPoints.append(RealPt)
            print('*********************************These are the SlbPts2:',SlbPts2)
        display(SlabPoints)
        ax.scatter(*zip(*SlabPoints),c='red')
    #****************************************Defines the holes per Slab*************************************************************
        Holes=[]
        if bool(i.HasOpenings) != False:
            for j in i.HasOpenings:
                Test=ifcopenshell.util.placement.get_local_placement(j.RelatedOpeningElement.ObjectPlacement)
                print('**********************************************************************')
                print(Test)
                print('+++++++++++++++++++++++++++++++++++++++++++++++++',j.RelatedOpeningElement.is_a('ifcOpeningElement'))
                if j.RelatedOpeningElement.Name == None:
                    try: 
                        #Convert the angle of the opening to the normalized angle of the host
                        try:
                            #directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios
                            directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios
                            
                        except:
                            print('No axis or refdirection found, using default placement')
                            directionX=(1,0,0)
                        Angle=-(angle_between((NormX[0],NormX[1]),(directionX[0],directionX[1])))
                        Coord=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.Location.Coordinates
                        PPoint=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.Location.Coordinates
                        print('these are the coordinates from the area of the opening: {} and these are the ones from the opening element itself {}'.format(Coord,PPoint))
                        Origin=(Coord[0]-PPoint[0],Coord[1]+PPoint[1])
                        Opn=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea
                        Pts=(-(Opn.XDim)/2,-(Opn.YDim)/2),(-(Opn.XDim)/2,(Opn.YDim)/2),((Opn.XDim)/2,(Opn.YDim)/2),((Opn.XDim)/2,-(Opn.YDim)/2)
                        NewPts=[] 
                        for k in Pts:
                            if j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]<0:
                                print('this is k:',k)
                                Ri=rotate_pt(k,Angle)
                                print('this is k after rotation, ',Ri)
                                NewPt=(Ri[0]-Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                print(PlacePoint[0],PlacePoint[1])
                                print('And this is the final point',NewPt)
                            elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]>=0:
                                print('this is k:',k)
                                Ri=rotate_pt(k,Angle)
                                print('this is k after rotation, ',Ri)
                                NewPt=(Ri[0]-Origin[0]+SlabLoc[0],Ri[1]-Origin[1]+SlabLoc[1])
                                print(PlacePoint[0],PlacePoint[1])
                                print('And this is the final point',NewPt)
                            NewPts.append(NewPt)
                        Holes.append(NewPts)
                        
                        print('These are the points of the void')
                        print(Pts)
                        print('These are thetransformed points {}',format(NewPts))
                        print('and this is the origin of the shape: {}'.format(PPoint))
                    except:
                        try:
                            try:
                                directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios
                            except:
                                print('No axis or refdirection found, using default placement')
                                directionX=(1,0,0)
                            Angle=-(angle_between((NormX[0],NormX[1]),(directionX[0],directionX[1])))
                            Coord=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.Location.Coordinates
                            PPoint=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.Location.Coordinates
                            print('these are the coordinates from the area of the opening: {} and these are the ones from the opening element itself {}'.format(Coord,PPoint))
                            Origin=(Coord[0]-PPoint[0],Coord[1]+PPoint[1])
                            Opn=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea
                            Pts=Opn.OuterCurve.Points
                            NewPts=[] 
                            for k in Pts:
                                if j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]<0:
                                    print('this is k:',k)
                                    Ri=rotate_pt(k,Angle)
                                    print('this is k after rotation, ',Ri)
                                    NewPt=(Ri[0]-Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                    print(PlacePoint[0],PlacePoint[1])
                                    print('And this is the final point',NewPt)
                                elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]>=0:
                                    print('this is k:',k)
                                    Ri=rotate_pt(k,Angle)
                                    print('this is k after rotation, ',Ri)
                                    NewPt=(Ri[0]-Origin[0]+SlabLoc[0],Ri[1]-Origin[1]+SlabLoc[1])
                                    print(PlacePoint[0],PlacePoint[1])
                                    print('And this is the final point',NewPt)
                                NewPts.append(NewPt)
                            Holes.append(NewPts)
                            print('These are the points of the void')
                            print(Pts)
                            print('These are thetransformed points {}',format(NewPts))
                            print('and this is the origin of the shape: {}'.format(PPoint))
                        except:
                            continue
                else:
                    try:
                        try:
                            directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios
                            #directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios
                            
                        except:
                            print('No axis or refdirection found, using default placement')
                            directionX=(1,0,0)
                        Angle=-(angle_between((NormX[0],NormX[1]),(directionX[0],directionX[1])))
                        Coord=j.RelatedOpeningElement.ObjectPlacement.RelativePlacement.Location.Coordinates

                        Opn=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea
                        Pts=(-(Opn.XDim)/2,-(Opn.YDim)/2),(-(Opn.XDim)/2,(Opn.YDim)/2),((Opn.XDim)/2,(Opn.YDim)/2),((Opn.XDim)/2,-(Opn.YDim)/2)
                        NewPts=[] 
                        for k in Pts:
                            if j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]<0:
                                print('this is k:',k)
                                Ri=rotate_pt(k,Angle)
                                print('this is k after rotation, ',Ri)
                                NewPt=(Ri[0]+Coord[0]+SlabLoc[0],Ri[1]+Coord[1]+SlabLoc[1])

                                print('And this is the final point',NewPt)
                            elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea.Position.RefDirection.DirectionRatios[1]>=0:
                                print('this is k:',k)
                                Ri=rotate_pt(k,Angle)
                                print('this is k after rotation, ',Ri)
                                NewPt=(Ri[0]+Coord[0]+SlabLoc[0],Ri[1]+Coord[1]+SlabLoc[1])

                                print('And this is the final point',NewPt)
                            NewPts.append(NewPt)
                        Holes.append(NewPts)

                        print('These are the points of the void')
                        print(Pts)
                        print('These are thetransformed points {}',format(NewPts))
                        print('and this is the origin of the shape: {}'.format(PPoint))
                    except Exception as e:
                        print('entro en la excepción deseada********************************************************************************')
                        print(e)
                        try:
                            try:
                                directionX=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios
                                
                            except Exception as e:
                                print(e,'No axis or refdirection found, using default placement')
                                directionX=(1,0,0)
                            Angle=(angle_between((NormX[0],NormX[1]),(directionX[0],directionX[1])))
                            #HostC=j.RelatedOpeningElement.ObjectPlacement.PlacementRelTo.PlacesObject[0].ObjectPlacement.RelativePlacement.
                            Coord=j.RelatedOpeningElement.ObjectPlacement.RelativePlacement.Location.Coordinates
                            PPoint=j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.Location.Coordinates
                            print('these are the coordinates from the area of the opening: {} and these are the ones from the opening element itself {}'.format(Coord,PPoint))
                            Origin=(Coord[0]+PPoint[0],Coord[1]+PPoint[1])
                            print('This is the Origin')
                            print(Origin)
                            Opn=j.RelatedOpeningElement.Representation.Representations[0].Items[0].SweptArea
                            Pts=Opn.OuterCurve.Points
                            print('These are the polycurve points',Pts)
                            NewPts=[] 
                            for k in Pts:
                                if j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.Axis.DirectionRatios[2]<0:
                                    if j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios[1]<0:
                                        print('this is k:',k)
                                        Pt=k.Coordinates
                                        Pt2=(Pt[0],-Pt[1])
                                        Ri=rotate_pt(Pt2,Angle)
                                        print('this is k after rotation, ',Ri)
                                        NewPt=(Ri[0]+Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                        print('this is the origin',Origin)
                                        print('And this is the final point',NewPt)
                                    elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios[1]>=0:
                                        print('this is k:',k)
                                        Pt=k.Coordinates
                                        Pt2=(Pt[0],-Pt[1])
                                        Ri=rotate_pt(Pt2,-Angle)
                                        print('this is k after rotation, ',Ri)
                                        NewPt=(Ri[0]+Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                        print('this is the origin',Origin)
                                        print('And this is the final point',NewPt)
                                    NewPts.append(NewPt)
                                elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.Axis.DirectionRatios[2]>=0:
                                    if j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios[1]<0:
                                        print('this is k:',k)
                                        Pt=k.Coordinates
                                        Pt2=(Pt[0],Pt[1])
                                        Ri=rotate_pt(Pt2,Angle)
                                        print('this is k after rotation, ',Ri)
                                        NewPt=(Ri[0]+Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                        print('this is the origin',Origin)
                                        print('And this is the final point',NewPt)
                                    elif j.RelatedOpeningElement.Representation.Representations[0].Items[0].Position.RefDirection.DirectionRatios[1]>=0:
                                        print('this is k:',k)
                                        Pt=k.Coordinates
                                        Pt2=(Pt[0],Pt[1])
                                        Ri=rotate_pt(Pt2,-Angle)
                                        print('this is k after rotation, ',Ri)
                                        NewPt=(Ri[0]+Origin[0]+SlabLoc[0],Ri[1]+Origin[1]+SlabLoc[1])
                                        print('this is the origin',Origin)
                                        print('And this is the final point',NewPt)
                                    NewPts.append(NewPt)
                            Holes.append(NewPts)
                            print('These are the points of the void')
                            print(Pts)
                            print('These are thetransformed points {}',format(NewPts))
                            print('and this is the origin of the shape: {}'.format(PPoint))
                        except Exception as e:
                            print(e)
                            continue
    #**************************************Define the representation of the holes in the slab***************************************************************        
            ax.scatter(*zip(*NewPts),c='blue')
            ptn=0
            for m in enumerate(NewPts):
                ax.annotate(ptn,(m[1][0],m[1][1]))
                ptn=ptn+1
        else:
            pass
    #**************************************Define the shapely shape of the slab trying to include the holes if they exist***********************************
        try:
            SlabShp= shapely.geometry.Polygon(SlabPoints,Holes)
        except:
            SlabShp= shapely.geometry.Polygon(SlabPoints)
        ax.add_patch(PolygonPatch(SlabShp, fc='red', alpha=0.2))
        try:
            print('Entro en el intento')
            SlabHoles=[]
            if len(Holes)>0:
                for n in Holes:
                    print(n)
                    VoidBound=shapely.geometry.Polygon(n)
                    SlabHoles.append(VoidBound)
            else:
                VoidBound=None
        except Exception as e:
            print(e)
            print('entro en la excepción')
            VoidBound=0
            pass
        print(SlabShp.area)
    #*******************************Defines the Slab entity*****************************************************************************************
        CurrSlab=Slab(i.GlobalId,SlabShp.area,SlabPoints,i.Tag,SlabShp,Holes,VoidBound)
        CurrSlab.Voids=SlabHoles
        CurrSlab.VoidPts=Holes
        SlabEn.append(CurrSlab)         
        Shapes.append(SlabShp)

    #********************************Display the shapes in the current set (level or total)*********************************************************
    plt.show()
    #***********************************************************************************************************************************************
    display(SlabEn)
    display(Shapes)
    #*********************************for each shape in the set, merge with the previous one********************************************************
    TempShp=None
    VoidShp=None
    for i in SlabEn:
        print('This is the Slab: {}, with Guid {} and a total Area of {}, the boundary points are: {} and the tag in authoring software is {}'.format(i,i.id,i.Area,i.BPoints,i.Tag))
        SlabShape=i.Shp
        VoidShape=i.Voids
        display(VoidShape)
        if TempShp==None:
            TempShp=SlabShape
            try:
                VoidShp=VoidShape[0]
            except:
                pass
        else:    
            TempShp=TempShp.union(SlabShape)
            if VoidShape==None:
                pass
            else:
                for j in VoidShape:
                    print('This is j:',j)
                    if VoidShp == None:
                        VoidShp=j
                    VoidShp=VoidShp.union(j)
        #TempShp=SlabShape.union(SlabEn[(i.index)+1].Shp)
    print('THIS IS THE TOTAL AREA OF SLAB OF THE SET (LEVEL OR SELECTION)',TempShp.area)
    LvlSlabShp[Lev.LongName]=TempShp
    #**********************************This is Wall Module!!!!****************************************************************************************
    #**********************************Find all the walls in the file and declare basic variables*****************************************************
    Walls= IfcFile.by_type('ifcWall')

    AllPoints=[]
    Vectors=[]
    Angles=[]
    WallAxPt=[]
    ConWalls=[]
    Location=[]
    PtLst=[]
    WallSet=[]
    WallEval=[]
    WallLoop=[]
    ItCount=0

    ExternalWalls=[]
    InternalWalls=[]
    for i in Walls:
        PsetCom=ifcopenshell.util.element.get_psets(i)['Pset_WallCommon']
        if PsetCom['IsExternal'] == True:
            ExternalWalls.append(i)
        else:
            InternalWalls.append(i)


    WallLev={}
    for i in ExternalWalls:
        if i.ContainedInStructure[0].RelatingStructure.LongName not in WallLev:
            WallLev[i.ContainedInStructure[0].RelatingStructure.LongName]=list()
            WallLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
        elif i.ContainedInStructure[0].RelatingStructure.LongName in WallLev:
            WallLev[i.ContainedInStructure[0].RelatingStructure.LongName].append(i)
        else:
            continue
    if bool(Lev)!= False:
        try:
            WallSet=WallLev[Lev.LongName]
        except:
            pass
    else:
        WallSet=ExternalWalls


    if bool(WallSet)!=False:
        for i in WallSet:
            WallEval.append(i)
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
        if bool(WallLoop) == False:
            continue
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
                        pass
                elif j is None:
                    OrdPts.append(i)
                elif i is None:
                    pass
                            
                    
                    
                
                        
            print('These are the ordinated pointsxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx')






            fig, ax = plt.subplots()

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
        #****************************************Creates the shape representation from the walls********************************************
        ax.scatter(*zip(*ExtOrd),c='green')
        ax.scatter(*zip(*TmpLst),c='cyan')

        ax.scatter(*zip(*IntOrd),c='black')
        print(Axs.wkt)
        print(Axs.is_valid)
        ptn=0
        for i in enumerate(IntOrd):
            ax.annotate(ptn,(i[1][0],i[1][1]))
            ptn=ptn+1
        ax.scatter(cent[0],cent[1],c='green')
        ax.add_patch(PolygonPatch(Axs, fc='red', alpha=0.2))
        ax.add_patch(PolygonPatch(ExtShp, fc='black', alpha=0.2))
        ax.add_patch(PolygonPatch(IntShp, fc='green', alpha=0.2))
        print(Axs.area)
        print('this is the floor area',ExtShp.area)
        plt.show()
  
    #********************************Combining all the shapes in the same graph**************************************************
        fig, ax = plt.subplots()
        for i in Shapes:
            ax.add_patch(PolygonPatch(i, fc='red', alpha=0.2))
        ax.add_patch(PolygonPatch(Axs, fc='red', alpha=0.2))
        for i in ExtShapes:
            ax.add_patch(PolygonPatch(i, fc='black', alpha=0.2))
        for i in IntShapes:
            ax.add_patch(PolygonPatch(i, fc='green', alpha=0.2))
        for i in PlotShps:
            ax.add_patch(PolygonPatch(i, fc='green', alpha=0.2))
        plt.show()
    #**********************************Merge the Area of the Walls with the Area of floor to get the overall shape******************************
    #fig, ax = plt.subplots()   
        for i in ExtShapes:
            TempShp=TempShp.union(i)
        '''try:
            TempShp=TempShp.difference(VoidShp) 
        except:
            pass'''
        LvlWallExtShp[Lev.LongName]=ExtShapes
        LvlWallIntShp[Lev.LongName]=IntShapes
        ax.add_patch(PolygonPatch(TempShp, fc='Blue', alpha=0.1))    
        plt.show()               
    #**********************************Creates a new representation of the merged shape************************************************************************
        fig, ax = plt.subplots()
        ax.add_patch(PolygonPatch(TempShp, fc='red', alpha=0.2))
        try:
            ax.add_patch(PolygonPatch(VoidShp, fc='blue', alpha=0.2))
        except:
            pass
        print('THIS IS THE TOTAL AREA OF THE SET (LEVEL {} OR SELECTION)'.format(Lev.LongName),TempShp.area)
        plt.show()
        AreaBrutaC+=TempShp.area
        LevABC[Lev.LongName]=TempShp.area
        LvlImpShp[Lev.LongName]=TempShp
        ImpShp.append(TempShp)
    else:
        fig, ax = plt.subplots()
        ax.add_patch(PolygonPatch(TempShp, fc='red', alpha=0.2))
        try:
            ax.add_patch(PolygonPatch(VoidShp, fc='blue', alpha=0.2))
        except:
            pass
        print('THIS IS THE TOTAL AREA OF THE SET (LEVEL {} OR SELECTION)'.format(Lev.LongName),TempShp.area)
        plt.show()
        AreaBrutaC+=TempShp.area
        LevABC[Lev.LongName]=TempShp.area
        LvlImpShp[Lev.LongName]=TempShp
        ImpShp.append(TempShp)
    #*****************************These are the uselful variables****************************
    '''print('These are the boundary walls and its points:')
    for i in WallObjects:
        print('This is the GUID')
        print(i.id)
        print('These are the wall Points')
        print(i.BPoints)
    print('These are the ordinated exterior points:')
    print(ExtOrd)
    print('These are the ordinated interior points:')
    print(IntOrd)
    print('These are the points of the site:')
    print(SitePoints[0])
    print('These are the Slabs used for the shape:')
    for i in SlabEn:
        print(i.id)
        print(i.BPoints)
    print('And these are the points of the merged shape:')
    print('These are the shape points*************************************************')
    #print(TempShp.exterior.coords.xy)
    print(list(zip(*TempShp.exterior.coords.xy)))
    #And these are the voids of the merged shape in case of need
    try:
        print('These are the void points**********************************')
        print(VoidShp.interior.coords.xy)
    except:
        pass'''
if Lev.LongName not in LvlImpShp:
    LvlImpShp[Lev.LongName]=TempShp
    ImpShp.append(TempShp)
fig, ax = plt.subplots()
ExclLvl=[]
exclusion='Y'
while exclusion != 'N':
    exclusion=input('Do you want to exclude any level shape from the implantation area?(Y/N)')
    while exclusion != 'N':
        print('Select the level to exclude:')
        for i in ChsnLvl:
            print(i.LongName)    
        ExclLvl.append(input(''))
        cont=input('Do you want to exclude another level?(Y/N)')
        if cont == 'N':
            exclusion = 'N'
display(ExclLvl)
ImpS=()
for i,j in LvlImpShp.items():
    if i in ExclLvl:
        continue
    elif bool(ImpS)==False:
        ImpS=j
    else:
        ImpS=ImpS.union(j)
ax.add_patch(PolygonPatch(ImpS, fc='Red', alpha=0.5))
plt.show()  
fig, ax = plt.subplots()     
Impz=()
for i in ImpShp:
    if bool(Impz)==False:
        Impz=i
    else:
        Impz=Impz.union(i)
ax.add_patch(PolygonPatch(Impz, fc='blue', alpha=0.5))
ax.add_patch(PolygonPatch(ImpS, fc='Red', alpha=0.1))


ExclLvl=[]
exclusion='Y'
while exclusion != 'N':
    exclusion=input('Do you want to exclude any level shape from the ABC area?(Y/N)')
    while exclusion != 'N':
        print('Select the level to exclude:')
        for i in ChsnLvl:
            print(i.LongName)    
        ExclLvl.append(input(''))
        cont=input('Do you want to exclude another level?(Y/N)')
        if cont == 'N':
            exclusion = 'N'
f = open("AreaImplicitExtractionResults.txt", 'w')
sys.stdout = f
AreaBrutaC=()
for i,j in LvlImpShp.items():
    if i in ExclLvl:
        continue
    elif bool(AreaBrutaC)==False:
        AreaBrutaC=j.area
        print('for the level',i,'the ABC area is',j.area)
    else:
        AreaBrutaC=AreaBrutaC+j.area
        print('for the level',i,'the ABC area is',j.area)
print('----------------------------------------------RESULTS--------------------------------------------------')
print('The plot area is:',AreaDoLote)
AreaDeImpl=ImpS.area
print('área de implantação is {}'.format(ImpS.area))
AreaDeLogra=AreaDoLote-AreaDeImpl
print('ÁreaDeLogradouro is',AreaDeLogra)
print('área de impermeabilização is {}'.format(Impz.area))
print('The Área Bruta is',AreaBrutaC)
f.close()
plt.show()
