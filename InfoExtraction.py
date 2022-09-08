import sys
from tkinter import Y
import ifcopenshell 
import ifcopenshell.util
import ifcopenshell.util.element
import ifcopenshell.util.selector
import sys
import ifcopenshell.geom

'''import ifcopenshell.geom.occ_utils
settings = ifcopenshell.geom.settings()
occ_display = ifcopenshell.geom.occ_utils.initialize_display()
#settings.set(settings.USE_PYTHON_OPENCASCADE, True)'''
def display(arg):
    for i in arg:
        print(i)

def Area(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area



file = ifcopenshell.open('C:\Académicos BIMAPLUS\BIM A+7\IfcOpenShell\TESTES\LV_EdificioHabitacional_R10.ifc')
Building=file.by_type('IfcBuilding')
spaces = file.by_type('IfcSpace')

ProjectInfo={}
BPsets=(ifcopenshell.util.element.get_psets(Building[0]))
try:
    DigBP=BPsets['BIMAPLUS_DBP']
    for i,j in DigBP.items():
        print(i,j)
        ProjectInfo[i]=j
except:
    print('No BIMAPLUS_DBP property set found')



print("-----------------------",spaces[0].get_info(),spaces[0].IsDecomposedBy)
print('-------------------------------------',len(spaces), "spaces in the project:---------------------------------------------")
print(spaces)


SpaceArea={}
InstCount=[]
AreaImp={}
ttlImp=float()
print('*****************************************Started looking for Area de implantacao***********************************************************')
for i in spaces:
    psets= (ifcopenshell.util.element.get_psets(i))
    try:
        DigBP=psets['BIMAPLUS_DBP']
    except:
        print('No BIMAPLUS_DBP property set found')
        continue
    if bool(DigBP['AreaType']) == True:
        try:        
            qto= psets['BaseQuantities']
            print ('Area found in:',i.LongName,'Area:',qto['GrossFloorArea'])
            if DigBP['AreaType'] not in SpaceArea:
                SpaceArea[DigBP['AreaType']]=qto['GrossFloorArea']
            elif DigBP['AreaType'] in SpaceArea:
                SpaceArea[DigBP['AreaType']]=SpaceArea[DigBP['AreaType']]+qto['GrossFloorArea']
            ttlImp=float(qto['GrossFloorArea'])+ttlImp
            continue
        except:
            print('No area was found in the space ',i.LongName,' properties')
            continue
print('********************************************this is the end of the execution***************************************************** Results: ')
f = open("AreaExtractionResults.txt", 'w')
sys.stdout = f
for i,j in ProjectInfo.items():
    print('This is the information attribute "{}", for this project it is:{}'.format(i,j))
print('----------------------------------------------------------------------------------------------------------------------------------------------------')
for i,j in SpaceArea.items():
    print('This is the area found as "{}" explicitly in the model:{} m2'.format(i,j))

f.close()
