from cgi import test
from logging import exception
from msilib.schema import Class
from operator import index
import os
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




NormX=(1,0,0)
Location=[]
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\SiteTest17.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\IFC_tt7.ifc')
#IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\COL_4BEDROOM DUPLEX.ifc')
IfcFile = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\LV_EdificioHabitacional_R07.ifc')

Levels=IfcFile.by_type('IfcBuildingStorey')
print('These are the levels of the project:')
for i in Levels:
    print(i.LongName)
SelLev=input('Please select a Level from the list above')

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
if bool(SelLev)!= False:
    SlabSet=SlabLev[SelLev]
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
print(len(SlabSet))
print(len(SlabEn))
#**********************************Creates a new representation of the merged shape************************************************************************
fig, ax = plt.subplots()
ax.add_patch(PolygonPatch(TempShp, fc='red', alpha=0.2))
try:
    ax.add_patch(PolygonPatch(VoidShp, fc='blue', alpha=0.2))
except:
    pass
print('THIS IS THE TOTAL AREA OF THE SET (LEVEL OR SELECTION)',TempShp.area)
plt.show()

