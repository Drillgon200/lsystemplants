import omni.ext
import omni.ui as ui
import omni.ui._ui
import carb.input
import omni.kit.mesh.raycast
import omni.graph.ui
import omni.kit.viewport_legacy
import omni.usd
from pxr import CameraUtil, Gf, UsdGeom, Vt, Sdf, UsdShade, Usd
from pxr.Gf import Frustum
import carb.dictionary
import carb.settings
from typing import List
from typing import Tuple
import random
import math
import os
import sys
import traceback

VIEWPORT_APERTURE_COMFORM_POLICY = "/app/hydra/aperture/conform"
LABEL_WIDTH = 120
VSPACING = 4
HSPACING = 12

# Any class derived from `omni.ext.IExt` in top level module (defined in `python.modules` of `extension.toml`) will be
# instantiated when extension gets enabled and `on_startup(ext_id)` will be called. Later when extension gets disabled
# on_shutdown() is called.
class MyExtension(omni.ext.IExt):

    class PlantSettings:
        
        def reset_general(self):
            self.branchSeed.as_int = -1
            self.leafSeed.as_int = -1

            self.endScale.as_float = 1
            self.width.as_float = 20
            self.length.as_float = 10
            self.maxIterations.as_int = 100
            self.upwardsPull.as_float = 0.4
            self.gravityPower.as_float = 1.0
            self.upDirection = Gf.Vec3d(0, 1, 0)
            self.branchLeafIterations.as_int = 3
            self.branchWidthScale.as_float = 0.25
            self.branchSizeDecay.as_float = 0.7
            self.lengthDecay.as_float = 0.9
        
        def reset_branch(self):
            self.subbranchChance.as_float = 0.3
            self.subbranchBranchChanceDecay.as_float = 0.5
            self.subbranchScaleMin.as_float = 5
            self.subbranchIterationMin.as_int = 3
            self.subbranchDownscale.as_float = 0.8
            self.subbranchAngleDefault.as_float = 55
            self.subbranchAngleMin.as_float = -5
            self.subbranchAngleMax.as_float = 15
            self.newBranchDirectionRange.as_float = 30
            self.branchOffValidRange.as_float = 140
            self.supportDecay.as_float = 0.8

        def reset_leaf(self):
            self.leafScale.as_float = 30
            self.leafDownturn.as_float = 0.6
            self.branchLeafScaleMax.as_float = 6
            self.leafChance.as_float = 0.2
            self.leafAngleDefault.as_float = 55
            self.leafAngleMin.as_float = -5
            self.leafAngleMax.as_float = 15

        def reset_textures(self):
            self.leafTexture.as_string = "leaf.png"
            self.branchTexture.as_string = "stem.png"

        def reset(self):
            self.reset_general()
            self.reset_branch()
            self.reset_leaf()
            self.reset_textures()

        def __init__(self):
            self.branchSeed = ui.SimpleIntModel()
            self.leafSeed = ui.SimpleIntModel()

            self.endScale = ui.SimpleFloatModel()
            self.width = ui.SimpleFloatModel()
            self.length = ui.SimpleFloatModel()
            self.maxIterations = ui.SimpleIntModel()
            self.upwardsPull = ui.SimpleFloatModel()
            self.gravityPower = ui.SimpleFloatModel()
            self.upDirection = Gf.Vec3d(0, 1, 0)
            self.branchLeafIterations = ui.SimpleIntModel()
            self.branchWidthScale = ui.SimpleFloatModel()
            self.branchSizeDecay = ui.SimpleFloatModel()
            self.lengthDecay = ui.SimpleFloatModel()

            self.subbranchChance = ui.SimpleFloatModel()
            self.subbranchBranchChanceDecay = ui.SimpleFloatModel()
            self.subbranchScaleMin = ui.SimpleFloatModel()
            self.subbranchIterationMin = ui.SimpleIntModel()
            self.subbranchDownscale = ui.SimpleFloatModel()
            self.subbranchAngleDefault = ui.SimpleFloatModel()
            self.subbranchAngleMin = ui.SimpleFloatModel()
            self.subbranchAngleMax = ui.SimpleFloatModel()
            self.newBranchDirectionRange = ui.SimpleFloatModel()
            self.branchOffValidRange = ui.SimpleFloatModel()
            self.supportDecay = ui.SimpleFloatModel()

            self.leafScale = ui.SimpleFloatModel()
            self.leafDownturn = ui.SimpleFloatModel()
            self.branchLeafScaleMax = ui.SimpleFloatModel()
            self.leafChance = ui.SimpleFloatModel()
            self.leafAngleDefault = ui.SimpleFloatModel()
            self.leafAngleMin = ui.SimpleFloatModel()
            self.leafAngleMax = ui.SimpleFloatModel()

            self.leafTexture = ui.SimpleStringModel()
            self.branchTexture = ui.SimpleStringModel()

            self.group = ui.SimpleStringModel("")

            self.reset()
        
        def serialize_to_prim(self, prim, origin, direction, branchSeed, leafSeed):
            prim.SetCustomDataByKey("PlantBSeed", branchSeed)
            prim.SetCustomDataByKey("PlantLSeed", leafSeed)

            prim.SetCustomDataByKey("PlantOX", origin[0])
            prim.SetCustomDataByKey("PlantOY", origin[1])
            prim.SetCustomDataByKey("PlantOZ", origin[2])
            prim.SetCustomDataByKey("PlantDX", direction[0])
            prim.SetCustomDataByKey("PlantDY", direction[1])
            prim.SetCustomDataByKey("PlantDZ", direction[2])

        def deserialize_from_prim(self, prim):
            branchSeed = prim.GetCustomDataByKey("PlantBSeed")
            leafSeed = prim.GetCustomDataByKey("PlantLSeed")

            ox = prim.GetCustomDataByKey("PlantOX")
            oy = prim.GetCustomDataByKey("PlantOY")
            oz = prim.GetCustomDataByKey("PlantOZ")
            dx = prim.GetCustomDataByKey("PlantDX")
            dy = prim.GetCustomDataByKey("PlantDY")
            dz = prim.GetCustomDataByKey("PlantDZ")

            return Gf.Vec3d(ox, oy, oz), Gf.Vec3d(dx, dy, dz), branchSeed, leafSeed


    def __init__(self):
        self.isDrawingLine = False
        self.isDrawingPlant = False
        self.linePointAnchor = None
        self.cableUseCurves = ui.SimpleBoolModel(True)
        self.cableSize = ui.SimpleFloatModel(40)
        self.plantSettings = self.PlantSettings()

    class Vertex:
        def __init__(self, pos: Gf.Vec3f, norm: Gf.Vec3f, tex: Gf.Vec2f):
            self.pos = pos
            self.norm = norm
            self.tex = tex

    def draw_line(self, pointA: Gf.Vec3d, pointB: Gf.Vec3d, firstWidth: float, secondWidth: float):
        stage = omni.usd.get_context().get_stage()

        curve = omni.usd.UsdGeom.BasisCurves.Define(stage, omni.usd.get_stage_next_free_path(stage, "/lplants/branch", True))
        curve.CreateCurveVertexCountsAttr([4])
        curve.CreatePointsAttr([pointA, pointA, pointB, pointB])
        curve.CreateWidthsAttr([firstWidth, firstWidth, secondWidth, secondWidth])
    
    def make_mat3d(self, a: Gf.Vec3d, b: Gf.Vec3d, c: Gf.Vec3d):
        return Gf.Matrix3d(a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2])
       # return Gf.Matrix3d(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2])

    #https://box2d.org/posts/2014/02/computing-a-basis/
    def compute_basis(self, a: Gf.Vec3d):
        if(abs(a[0]) >= 0.57735):
            b = Gf.Vec3d(a[1], -a[0], 0)
        else:
            b = Gf.Vec3d(0, a[2], -a[1])
        b = b.GetNormalized()
        return self.make_mat3d(b, a, Gf.Cross(a, b))

    def eval_bezier(self, pointA: Gf.Vec3d, ctrlA: Gf.Vec3d, ctrlB: Gf.Vec3d, pointB: Gf.Vec3d, t: float):
        tA = (ctrlA - pointA) * t + pointA
        tB = (ctrlB - ctrlA) * t + ctrlA
        tC = (pointB - ctrlB) * t + ctrlB
        tD = (tB - tA) * t + tA
        tE = (tC - tB) * t + tB
        tangent = tE - tD
        point = tangent * t + tD
        return point, tangent.GetNormalized()

    def eval_bezier_scalar(self, a: float, b: float, c: float, d: float, t: float):
        tA = (b - a) * t + a
        tB = (c - b) * t + b
        tC = (d - c) * t + c
        tD = (tB - tA) * t + tA
        tE = (tC - tB) * t + tB
        return (tE - tD) * t + tD
    
    def approximate_bezier_length(self, pointA: Gf.Vec3d, ctrlA: Gf.Vec3d, ctrlB: Gf.Vec3d, pointB: Gf.Vec3d):
        return Gf.GetLength(pointA - ctrlA) + Gf.GetLength(ctrlA - ctrlB) + Gf.GetLength(ctrlB - pointB)

    def draw_bezier(self, points: List[Gf.Vec3d], sizes: List[float], minV: float, maxV: float, geo: List[Vertex]):
        if(len(points) < 4):
            return
        #basis = self.compute_basis(Gf.Vec3d(points[1] - points[0]).GetNormalized())
        basis = Gf.Matrix3d(Gf.Rotation(Gf.Vec3d(0, 1, 0), Gf.Vec3d(points[1] - points[0]).GetNormalized()))
        previousCircle = []
        previousPoint = Gf.Vec3d(0, 0, 0)
        circleSegments = min(12, max(3, int(sizes[0] * 2)))
        geoPos = len(geo)
        totalLength = 0
        for i in range((len(points)-1)//3):
            pointA = points[i * 3]
            ctrlA = points[i * 3 + 1]
            ctrlB = points[i * 3 + 2]
            pointB = points[i * 3 + 3]
            approxLength = self.approximate_bezier_length(pointA, ctrlA, ctrlB, pointB)
            totalLength += approxLength
        currentV = minV
        dV = maxV - minV
        vScale = totalLength / (2 * math.pi * (sizes[0] * 0.5))
        for i in range((len(points)-1)//3):
            previousCircle.clear()

            pointA = points[i * 3]
            ctrlA = points[i * 3 + 1]
            ctrlB = points[i * 3 + 2]
            pointB = points[i * 3 + 3]
            sizeA = sizes[i * 3]
            sizeB = sizes[i * 3 + 1]
            sizeC = sizes[i * 3 + 2]
            sizeD = sizes[i * 3 + 3]

            approxLength = self.approximate_bezier_length(pointA, ctrlA, ctrlB, pointB)
            substeps = int(approxLength)
            substeps = min(20, max(2, substeps))
            
            nextV = currentV + (approxLength/totalLength) * dV
            prevV = currentV

            for i in range(substeps):
                t = i / (substeps-1)
                point, tangent = self.eval_bezier(pointA, ctrlA, ctrlB, pointB, t)
                point = Gf.Vec3d(point)
                tangent = Gf.Vec3d(tangent)

                v = ((nextV - currentV) * t + currentV) * vScale

                nextBasisVec = Gf.Cross(tangent, basis.GetColumn(0)).GetNormalized()
                basis = self.make_mat3d(Gf.Cross(nextBasisVec, tangent).GetNormalized(), tangent, nextBasisVec)

                size = self.eval_bezier_scalar(sizeA, sizeB, sizeC, sizeD, t)
                circle = []
                for c in range(circleSegments):
                    angle = (c/circleSegments) * 2 * math.pi
                    circlePos = Gf.Vec3d(math.cos(angle), 0, math.sin(angle)) * size * 0.5
                    circlePos = basis * circlePos + point
                    circle.append(circlePos)
                
                for c in range(len(previousCircle)):
                    # Generate geometry
                    cur = c
                    next = (c + 1) % circleSegments
                    curU = c / circleSegments
                    nextU = (c + 1) / circleSegments
                    geo.append(self.Vertex(Gf.Vec3f(previousCircle[cur]), Gf.Vec3f((previousCircle[cur] - previousPoint).GetNormalized()), Gf.Vec2f(curU, prevV)))
                    geo.append(self.Vertex(Gf.Vec3f(previousCircle[next]), Gf.Vec3f((previousCircle[next] - previousPoint).GetNormalized()), Gf.Vec2f(nextU, prevV)))
                    geo.append(self.Vertex(Gf.Vec3f(circle[next]), Gf.Vec3f((circle[next] - point).GetNormalized()), Gf.Vec2f(nextU, v)))
                    geo.append(self.Vertex(Gf.Vec3f(circle[cur]), Gf.Vec3f((circle[cur] - point).GetNormalized()), Gf.Vec2f(curU, v)))
                
                previousCircle = circle
                previousPoint = point
                prevV = v
            
            currentV = nextV
        
        # Calculate normals
        #sameVertices = {}
        #i = geoPos
        #while i < len(geo):
        #    vert = geo[i]
        #    list = sameVertices.get(vert.pos)
        #    if list is None:
        #        list = []
        #        sameVertices[vert.pos] = list
        #    list.append(vert)
        #    geoPos += 1

        #for vertList in sameVertices:
        #    norm = Gf.Vec3f(0, 0, 0)
        #    for vert in vertList:
        #        norm = norm + vertList.normal
        #    norm = norm.GetNormalized()
        #    for vert in vertList:
        #        vert.normal = norm


                


    def draw_line_strip(self, rootNormal: Gf.Vec3d, lines: List[Tuple[float, Gf.Vec3d]], pathBase: str, geo: List[Vertex]):
        if(len(lines) < 2):
            return

        curvepoints = []
        sizes = []

        bezierControl = 0.35
        first = lines[0][1]
        curvepoints.append(first)
        #curvepoints.append((lines[1][1] - first) * bezierControl + first)
        curvepoints.append(first + rootNormal * Gf.GetLength(lines[1][1] - first) * bezierControl)
        sizes.append(lines[0][0])
        sizes.append(Gf.Lerp(bezierControl, lines[1][0], lines[0][0]))

        i = 1
        while i < len(lines) - 2:
            last = lines[i-1][1]
            cur = lines[i][1]
            nex = lines[i + 1][1]
            curToLast = (last - cur).GetNormalized()
            curToNext = (nex - cur).GetNormalized()
            tangent = Gf.Cross(curToLast + curToNext, Gf.Cross(curToLast, curToNext)).GetNormalized()
            lastLen = Gf.GetLength(last - cur)
            nextLen = Gf.GetLength(nex - cur)
            curvepoints.append(cur + tangent * lastLen * bezierControl)
            curvepoints.append(cur)
            curvepoints.append(cur - tangent * nextLen * bezierControl)
            sizes.append(Gf.Lerp(bezierControl, lines[i-1][0], lines[i][0]))
            sizes.append(lines[i][0])
            sizes.append(Gf.Lerp(bezierControl, lines[i+1][0], lines[i][0]))
            i += 1

        last = lines[-1][1]
        curvepoints.append((lines[-2][1] - last) * bezierControl + last)
        curvepoints.append(last)
        sizes.append(Gf.Lerp(bezierControl, lines[-2][0], lines[-1][0]))
        sizes.append(lines[-1][0])

        # Curves don't really have good enough UV coords set up
        self.draw_bezier(curvepoints, sizes, 0, 1, geo)

        #stage = omni.usd.get_context().get_stage()
        #curve = omni.usd.UsdGeom.BasisCurves.Define(stage, omni.usd.get_stage_next_free_path(stage, pathBase + "/branch", True))
        #curve.CreateCurveVertexCountsAttr([len(curvepoints)])
        #curve.CreatePointsAttr(curvepoints)
        #curve.CreateWidthsAttr(sizes)

    def add_quad_geometry(self, prevMesh, geo: List[Vertex], pathBase: str, name: str, texture: str):
        positions = []
        normals = []
        texCoords = []
        indices = []
        faceVertexCounts = []
        vertIdx = -1
        indexBuilder = {}

        for vert in geo:
            idx = indexBuilder.get(vert.pos)
            if idx is not None:
                indices.append(idx)
            else:
                vertIdx += 1
                positions.append(vert.pos)
                indices.append(vertIdx)

                indexBuilder[vert] = vertIdx
            # Apparently indices only work for positions???
            # I'm used to the indices being for the whole vertex
            normals.append(vert.norm)
            texCoords.append(vert.tex)
            
        for _ in range(len(indices)//4):
            faceVertexCounts.append(4)

        
        
        if(prevMesh is not None):
            mesh = prevMesh
        else:
            stage = omni.usd.get_context().get_stage()
            mesh = omni.usd.UsdGeom.Mesh.Define(stage, omni.usd.get_stage_next_free_path(stage, f"{pathBase}/{name}", True))
        mesh.GetPointsAttr().Set(Vt.Vec3fArray(positions))
        mesh.GetNormalsAttr().Set(Vt.Vec3fArray(normals))
        mesh.GetFaceVertexIndicesAttr().Set(indices)
        mesh.GetFaceVertexCountsAttr().Set(faceVertexCounts)
        mesh.SetNormalsInterpolation("faceVarying")

        prim = mesh.GetPrim()
        
        sts_primvar = mesh.CreatePrimvar("st", Sdf.ValueTypeNames.Float2Array)
        sts_primvar.SetInterpolation("faceVarying")
        sts_primvar.Set(Vt.Vec2fArray(texCoords))
        mesh.CreateSubdivisionSchemeAttr("none")

        #omni.usd.get_context().get_selection().set_selected_prim_paths([prim.GetPath().pathString], False)

        
        
        if len(texture) != 0 and prevMesh is None:
            # Set material
            # Thanks Mati

            texPath = os.path.join(__file__, f"../resources/{texture}")
            test_mtl =f"lplant_mtl_{texture.split('.')[0]}"
            
            shaderPath = f"/World/Looks/{test_mtl}"
            if not omni.usd.is_path_valid(Sdf.Path(shaderPath)):
                # Create a Shader
                success, result = omni.kit.commands.execute('CreateAndBindMdlMaterialFromLibrary',mdl_name='OmniPBR.mdl',mtl_name='OmniPBR', prim_name=test_mtl, mtl_created_list=None,bind_selected_prims=False)

                # Get the Shader
                shader_prim = stage.GetPrimAtPath(shaderPath + "/Shader")
                shader_prim.CreateAttribute("inputs:diffuse_texture", Sdf.ValueTypeNames.Asset).Set(texPath)

                shader_prim.CreateAttribute("inputs:reflection_roughness_constant", Sdf.ValueTypeNames.Float).Set(1)
                shader_prim.CreateAttribute("inputs:specular_level", Sdf.ValueTypeNames.Float).Set(0)

                shader_prim.CreateAttribute("inputs:ao_to_diffuse", Sdf.ValueTypeNames.Float).Set(0.2)

                shader_prim.CreateAttribute("inputs:enable_opacity", Sdf.ValueTypeNames.Bool).Set(True)
                shader_prim.CreateAttribute("inputs:enable_opacity_texture", Sdf.ValueTypeNames.Bool).Set(True)
                shader_prim.CreateAttribute("inputs:opacity_texture", Sdf.ValueTypeNames.Asset).Set(texPath)
                shader_prim.CreateAttribute("inputs:opacity_mode", Sdf.ValueTypeNames.Int).Set(0) # mono_alpha
            omni.kit.commands.execute('BindMaterialCommand', prim_path=prim.GetPath().pathString, material_path=shaderPath, strength=UsdShade.Tokens.strongerThanDescendants) 
        return prim

    def add_leaf(self, options: PlantSettings, origin: Gf.Vec3d, direction: Gf.Vec3d, scale: float, rand: random.Random, geo: List[Vertex]):
        direction = (direction - options.upDirection * options.leafDownturn.as_float * (Gf.Clamp(Gf.Dot(direction, options.upDirection), 0, 1))).GetNormalized()
        leafAxis = Gf.Cross(options.upDirection, direction).GetNormalized()
        normal = Gf.Cross(direction, leafAxis).GetNormalized()
        lastNormal = normal
        vertexPair = [
            self.Vertex(Gf.Vec3f(origin + leafAxis * scale * 0.5), Gf.Vec3f(normal), Gf.Vec2f(0, 0)),
            self.Vertex(Gf.Vec3f(origin - leafAxis * scale * 0.5), Gf.Vec3f(normal), Gf.Vec2f(1, 0))
        ]
        geo.append(vertexPair[0])
        geo.append(vertexPair[1])
        rotations = [10, 15, 20]
        for i in range(3):
            rot = -rotations[i]
            rotationMat = Gf.Matrix3d(Gf.Rotation(leafAxis, rot))
            lastNormal = normal
            normal = (rotationMat * normal).GetNormalized()
            
            vertNormal = (lastNormal + normal).GetNormalized()
            
            v = (i + 1) / 4
            vertexPair = [
                self.Vertex(vertexPair[0].pos + Gf.Vec3f(direction * scale * 0.25), Gf.Vec3f(vertNormal), Gf.Vec2f(0, v)),
                self.Vertex(vertexPair[1].pos + Gf.Vec3f(direction * scale * 0.25), Gf.Vec3f(vertNormal), Gf.Vec2f(1, v))
            ]

            direction = (rotationMat * direction).GetNormalized()
            
            geo.append(vertexPair[1])
            geo.append(vertexPair[0])
            geo.append(vertexPair[0])
            geo.append(vertexPair[1])

        vertexPair = [
            self.Vertex(geo[-2].pos + Gf.Vec3f(direction * scale * 0.25), Gf.Vec3f(normal), Gf.Vec2f(0, 1)),
            self.Vertex(geo[-1].pos + Gf.Vec3f(direction * scale * 0.25), Gf.Vec3f(normal), Gf.Vec2f(1, 1))
        ]
        geo.append(vertexPair[1])
        geo.append(vertexPair[0])
        

    def generate_ls_plant(self, options: PlantSettings, origin: Gf.Vec3d, orientation: Gf.Matrix3d, scale: float, lengthScale: float, branchChance: float, support: float, iterations: int, branchRand: random.Random, leafRand: random.Random, linesList: List[List[Tuple[float, Gf.Vec3d]]], lines: List[Tuple[float, Gf.Vec3d]], leafGeo: List[Vertex]):
        direction: Gf.Vec3d = orientation.GetRow(1)
        if(iterations < 0 or scale < options.endScale.as_float):
            self.add_leaf(options, origin, direction, options.leafScale.as_float, leafRand, leafGeo)
            return
        
        newDirectionPitch = branchRand.uniform(-options.newBranchDirectionRange.as_float, options.newBranchDirectionRange.as_float)
        newDirection = Gf.Matrix3d(Gf.Rotation(orientation.GetRow(0), newDirectionPitch)) * direction
        newDirectionOrientation = branchRand.uniform(0, 360)
        newDirection = Gf.Matrix3d(Gf.Rotation(direction, newDirectionOrientation)) * newDirection
        newDirection += options.upDirection * options.upwardsPull.as_float #* (1-math.pow(Gf.Clamp(Gf.Dot(direction, up), 0, 1), 4))
        newDirection += -options.upDirection * (1 - support) * options.gravityPower.as_float
        newDirection = newDirection.GetNormalized()

        point2 = origin + newDirection * lengthScale

        lines.append((scale * options.branchWidthScale.as_float, point2))

        for _ in range(options.branchLeafIterations.as_int):
            branchRandom = branchRand.random()
            leafRandom = leafRand.random()
            if(branchRandom < branchChance and iterations > options.subbranchIterationMin.as_int and scale > options.subbranchScaleMin.as_float):
                randomOrigin = (point2 - origin) * branchRand.random() + origin
                #angleOldToNew = Gf.RadiansToDegrees(math.acos(Gf.Dot(direction, newDirection)))
                

                branchAngle = branchRand.uniform(options.subbranchAngleMin.as_float, options.subbranchAngleMax.as_float) + math.copysign(options.subbranchAngleDefault.as_float, -newDirectionPitch)
                orientationOpposite = (newDirectionOrientation + 180 + branchRand.uniform(-options.branchOffValidRange.as_float, options.branchOffValidRange.as_float)) % 360
                branchDirection = Gf.Matrix3d(Gf.Rotation(orientation.GetRow(0), branchAngle)) * direction
                newDirection = Gf.Matrix3d(Gf.Rotation(direction, orientationOpposite)) * newDirection

                branchList = []
                branchScale = scale * options.subbranchDownscale.as_float
                branchList.append((branchScale * options.branchWidthScale.as_float, randomOrigin))
                self.generate_ls_plant(options, randomOrigin, Gf.Matrix3d(Gf.Rotation(options.upDirection, branchDirection)), branchScale * options.branchSizeDecay.as_float, lengthScale * options.lengthDecay.as_float, branchChance * options.subbranchBranchChanceDecay.as_float, 1.0, iterations, branchRand, leafRand, linesList, branchList, leafGeo)
                linesList.append(branchList)
            elif(leafRandom < options.leafChance.as_float and scale < options.branchLeafScaleMax.as_float):
                randomOrigin = (point2 - origin) * leafRand.random() + origin
                branchAngle = leafRand.uniform(options.leafAngleMin.as_float, options.leafAngleMax.as_float) + math.copysign(options.leafAngleDefault.as_float, -newDirectionPitch)
                orientationOpposite = (newDirectionOrientation + 180 + leafRand.uniform(-options.branchOffValidRange.as_float, options.branchOffValidRange.as_float)) % 360
                branchDirection = Gf.Matrix3d(Gf.Rotation(orientation.GetRow(0), branchAngle)) * direction
                newDirection = Gf.Matrix3d(Gf.Rotation(direction, orientationOpposite)) * newDirection
                self.add_leaf(options, randomOrigin, branchDirection, options.leafScale.as_float, leafRand, leafGeo)
        
        self.generate_ls_plant(options, point2, Gf.Matrix3d(Gf.Rotation(options.upDirection, newDirection)), scale * options.branchSizeDecay.as_float, lengthScale, branchChance, support * options.supportDecay.as_float, iterations - 1, branchRand, leafRand, linesList, lines, leafGeo)

    def generate_plant_impl(self, options, branchMesh, leafMesh, origin, direction, branchSeed, leafSeed, pathBase):
        orientation = Gf.Matrix3d(Gf.Rotation(options.upDirection, direction))
        lines = [[]]
        lines[0].append((options.width.as_float * options.branchWidthScale.as_float, origin))
        leafGeo = []
        stemGeo = []
        
        rand = random.Random(branchSeed)
        branchRand = random.Random(rand.randint(0, sys.maxsize))
        rand.seed(leafSeed)
        leafRand = random.Random(rand.randint(0, sys.maxsize))
        
        self.generate_ls_plant(options, origin, orientation, options.width.as_float * 0.8, options.length.as_float, options.subbranchChance.as_float, 1.0, options.maxIterations.as_int, branchRand, leafRand, lines, lines[0], leafGeo)
        
        for lineSegment in lines:
            self.draw_line_strip(direction, lineSegment, pathBase, stemGeo)
        self.add_quad_geometry(leafMesh, leafGeo, pathBase, "leaves", options.leafTexture.as_string)
        stemPrim = self.add_quad_geometry(branchMesh, stemGeo, pathBase, "stem", options.branchTexture.as_string)
        if branchMesh is None:
            options.serialize_to_prim(stemPrim, origin, direction, branchSeed, leafSeed)

    def generate_plant(self, origin, direction, randSeed):
        options = self.plantSettings

        group = self.plantSettings.group.as_string
        group = ("/" + group) if group else ""
        stage = omni.usd.get_context().get_stage()
        pathBase = omni.usd.get_stage_next_free_path(stage, f"/lplants{group}/plant", True)

        branchSeed = options.leafSeed.as_int if options.leafSeed.as_int != -1 else (randSeed)
        leafSeed = options.leafSeed.as_int if options.leafSeed.as_int != -1 else (randSeed + 1)

        self.generate_plant_impl(options, None, None, origin, direction, branchSeed, leafSeed, pathBase)
    
    def regenerate_plant(self, stemMesh, leafMesh):
        options = self.plantSettings
        origin, direction, branchSeed, leafSeed = options.deserialize_from_prim(stemMesh.GetPrim())
        #print("Regenerating plant!")
        self.generate_plant_impl(options, stemMesh, leafMesh, origin, direction, branchSeed, leafSeed, None)

    def regenerate_selected_plants(self):
        stage = omni.usd.get_context().get_stage()
        selected = omni.usd.get_context().get_selection().get_selected_prim_paths()
        meshPaths = set()
        for path in selected:
            for prim in Usd.PrimRange(omni.usd.get_prim_at_path(Sdf.Path(path))):
                if prim.GetPrimTypeInfo().GetTypeName() == "Mesh":
                    meshPaths.add(str(prim.GetPrimPath()))
        for path in meshPaths:
            isPlantStem = path.endswith("stem") and "lplants" in path
            if not isPlantStem:
                continue
            leafPath = path[:path.rindex('/') + 1] + "leaves"
            stemPrim = omni.usd.get_prim_at_path(Sdf.Path(path))
            leafPrim = omni.usd.get_prim_at_path(Sdf.Path(leafPath))
            if not leafPrim.IsValid():
                continue
            if stemPrim.GetPrimTypeInfo().GetTypeName() == "Mesh" and leafPrim.GetPrimTypeInfo().GetTypeName() == "Mesh":
                stemMesh = UsdGeom.Mesh.Get(stage, Sdf.Path(path))
                leafMesh = UsdGeom.Mesh.Get(stage, Sdf.Path(leafPath))

                self.regenerate_plant(stemMesh, leafMesh)

    #This method in the brush widget class didn't do it right, so I guess I'm doing this myself
    def get_position_in_viewport(self, coords):
        if coords.x < 0.0 or coords.x > 1.0 or coords.y < 0.0 or coords.y > 1.0:
            return (None, None)

        x = coords.x * ui.Workspace.get_main_window_width()
        y = coords.y * ui.Workspace.get_main_window_height()

        windows = ui.Workspace.get_windows()
        #Ought to be more accurate than getting whatever window is named "Viewport"
        #There has to be a better way to get the topmost widow the mouse is over.
        #I do not know it, so iterating over all the windows it is!
        #There shouldn't be more than a few tens of windows at any given time anyway
        mouseOverWindow = None
        #I hope python's reversed function gives me an iterator and doesn't make a whole new list just to
        #iterate over it backwards
        for window in reversed(windows):
            #Check if the window is visible, the current one in dock, and the mouse is in it
            #Go with the first one because the reverse of the windows list is sorted front to back
            if window.visible and window.is_selected_in_dock() and x >= window.position_x and x <= (window.position_x + window.width) and y >= window.position_y and y <= (window.position_y + window.height):
                mouseOverWindow = window
                break
        
        #Do I really have to check the title here? That seems really unstable, but it's what other addons do
        if mouseOverWindow is None or "Viewport" not in mouseOverWindow.title:
            return (None, None)

        viewportBorderSize = mouseOverWindow.padding_x
        tab_bar_height = 22 if mouseOverWindow.dock_tab_bar_visible else 0

        windowXMin = mouseOverWindow.position_x + mouseOverWindow.padding_x
        windowXMax = mouseOverWindow.position_x + mouseOverWindow.width - mouseOverWindow.padding_x
        windowYMin = mouseOverWindow.position_y + mouseOverWindow.padding_y + tab_bar_height
        windowYMax = mouseOverWindow.position_y + mouseOverWindow.height - mouseOverWindow.padding_y

        if(x < windowXMin or x > windowXMax or y < windowYMin or y > windowYMax):
            return (None, None)
        x = min(max(x, windowXMin), windowXMax)
        y = min(max(y, windowYMin), windowYMax)

        x = (x - windowXMin) / (windowXMax - windowXMin)
        y = (y - windowYMin) / (windowYMax - windowYMin)

        return (x * 2 - 1, -(y * 2 - 1))

    def raycast_mouse(self, mX, mY):
        if mX is None or mY is None:
            return None
        viewport = omni.kit.viewport_legacy.get_viewport_interface()
        viewportWindow = viewport.get_viewport_window()
        active_camera = viewportWindow.get_active_camera()
        cam_prim = omni.usd.get_context().get_stage().GetPrimAtPath(active_camera)
        if not cam_prim.IsValid():
            return None

        gf_camera = UsdGeom.Camera(cam_prim).GetCamera()
        frustum = gf_camera.frustum
        resolution = viewportWindow.get_texture_resolution()
        
        CameraUtil.ConformWindow(frustum, CameraUtil.Crop, resolution[0] / resolution[1])

        ray = frustum.ComputePickRay(Gf.Vec2d(mX, mY))

        raycaster = omni.kit.mesh.raycast.get_mesh_raycast_interface()
        hitPoint = raycaster.closestRaycast(ray.startPoint, ray.direction, 100000)
        return hitPoint

    def mouse_event_handler(self, mouseInput):
        event = mouseInput.event
        
        if self.isDrawingLine:
            if event.type == carb.input.MouseEventType.LEFT_BUTTON_DOWN:
                x, y = self.get_position_in_viewport(event.normalized_coords)
                if(x is None):
                    self.isDrawingLine = False
                    self.linePointAnchor = None
                    return True
                point = self.raycast_mouse(x, y)

                if point is not None:
                    if self.linePointAnchor is not None:
                        pointA = Gf.Vec3f(self.linePointAnchor.position.x, self.linePointAnchor.position.y, self.linePointAnchor.position.z)
                        pointB = Gf.Vec3f(point.position.x, point.position.y, point.position.z)
                        normA = Gf.Vec3f(self.linePointAnchor.normal.x, self.linePointAnchor.normal.y, self.linePointAnchor.normal.z)
                        normB = Gf.Vec3f(point.normal.x, point.normal.y, point.normal.z)
                        line = pointB - pointA
                        lineLength = Gf.GetLength(line)
                        middleDown = (pointA + pointB) * 0.5 - Gf.Vec3f(0, lineLength * 0.2, 0)

                        points = [pointA, pointA + normA * 10.0, middleDown - line * 0.3, middleDown, middleDown + line * 0.3, pointB + normB * 10.0, pointB]
                        lineWidth = self.cableSize.as_float
                        widths = [lineWidth, lineWidth, lineWidth, lineWidth, lineWidth, lineWidth, lineWidth]

                        if not self.cableUseCurves.as_bool:
                            geo = []
                            self.draw_bezier(points, widths, 0, 1, geo)
                            self.add_quad_geometry(None, geo, "lplants/cables", "cable", "")
                        else:
                            stage = omni.usd.get_context().get_stage()
                            curve = omni.usd.UsdGeom.BasisCurves.Define(stage, omni.usd.get_stage_next_free_path(stage, "lplants/cables/cable", True))
                            curve.CreateCurveVertexCountsAttr([7])
                            curve.CreatePointsAttr(points)
                            curve.CreateWidthsAttr(widths)

                        self.linePointAnchor = None
                    else:
                        self.linePointAnchor = point
                return False
        elif self.isDrawingPlant:
            if event.type == carb.input.MouseEventType.LEFT_BUTTON_DOWN:
                x, y = self.get_position_in_viewport(event.normalized_coords)
                if(x is None):
                    self.isDrawingPlant = False
                    return True
                point = self.raycast_mouse(x, y)
                self.generate_plant(Gf.Vec3d(point.position.x, point.position.y, point.position.z), Gf.Vec3d(point.normal.x, point.normal.y, point.normal.z), random.randint(0, sys.maxsize))
                return False
        return True

    def key_event_handler(self, keyInput):
        if keyInput.event.type == carb.input.KeyboardEventType.KEY_PRESS:
            if(keyInput.event.input == carb.input.KeyboardInput.ESCAPE):
                self.isDrawingLine = False
                self.linePointAnchor = None
                self.isDrawingPlant = False
        return True

    # The mouse_moved_fn function takes 4 arguments.
    # I don't have a clue what they are, and I don't need them
    def on_plant_settings_modified(self, a, b, c, d):
        # WHAT IN THE WORLD IS d??? I JUST SPENT AN HOUR TRYING TO FIGURE OUT WHY
        # THIS FUNCTION WAS GETTING CALLED AT RANDOM TIMES (meaning when right mouse was released anywhere)
        # AND APPARENTLY d IS 0 WHEN THAT HAPPENS???
        # Oh well, guess I'll just do a check for that
        if d != 0:
            self.regenerate_selected_plants()

    def build_plant_general_settings(self):
        with ui.VStack(height=0, spacing=VSPACING):
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Seed", name="attribute_name", width=LABEL_WIDTH)
                ui.IntDrag(model=self.plantSettings.branchSeed, min=-10000, max=10000)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Seed", name="attribute_name", width=LABEL_WIDTH)
                ui.IntDrag(model=self.plantSettings.leafSeed, min=-10000, max=10000)
            with ui.HStack(spacing=HSPACING):
                ui.Label("End Scale", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.endScale, min=0.5, max=20, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Width", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.width, min=0.5, max=500, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Length", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.length, min=0.5, max=500, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Max Iterations", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.maxIterations, min=1, max=150, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Upwards Pull", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.upwardsPull, min=0, max=10, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Gravity Power", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.gravityPower, min=0, max=10, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Leaf Iterations", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.branchLeafIterations, min=0, max=15, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Width Scale", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.branchWidthScale, min=0, max=10, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Size Decay", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.branchSizeDecay, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Length Decay", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(self.plantSettings.lengthDecay, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            
            def reset_general_settings():
                self.plantSettings.reset_general()

            ui.Button("Reset General Settings", clicked_fn=lambda: reset_general_settings())
    
    def build_plant_branch_settings(self):
        with ui.VStack(height=0, spacing=VSPACING):
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Chance", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchChance, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Chance Decay", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchBranchChanceDecay, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Scale Min", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchScaleMin, min=0, max=500, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Iteration Min", name="attribute_name", width=LABEL_WIDTH)
                ui.IntDrag(model=self.plantSettings.subbranchIterationMin, min=0, max=150, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Downscale", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchDownscale, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Starting Angle", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchAngleDefault, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Angle Low Variation", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchAngleMin, min=-180, max=0, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branching Angle High Variation", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.subbranchAngleMax, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("New Branch Direction Range", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.newBranchDirectionRange, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Branch Off Valid Range", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.branchOffValidRange, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Support Decay", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.supportDecay, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)

            def reset_branch_settings():
                self.plantSettings.reset_branch()

            ui.Button("Reset Branch Settings", clicked_fn=lambda: reset_branch_settings())

    def build_plant_leaf_settings(self):
        with ui.VStack(height=0, spacing=VSPACING):
            
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Scale", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafScale, min=0, max=100, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Downturn", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafDownturn, min=0, max=10, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Max Branch Scale For Leaves", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.branchLeafScaleMax, min=0, max=500, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Chance", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafChance, min=0, max=1, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Starting Angle", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafAngleDefault, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf Low Variation", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafAngleMin, min=-180, max=0, mouse_released_fn=self.on_plant_settings_modified)
            with ui.HStack(spacing=HSPACING):
                ui.Label("Leaf High Variation", name="attribute_name", width=LABEL_WIDTH)
                ui.FloatSlider(model=self.plantSettings.leafAngleMax, min=0, max=180, mouse_released_fn=self.on_plant_settings_modified)

            def reset_leaf_settings():
                self.plantSettings.reset_leaf()
                
            ui.Button("Reset Leaf Settings", clicked_fn=lambda: reset_leaf_settings())

    def build_plant_settings(self):
        with ui.CollapsableFrame("General Settings", name="group"):
            self.build_plant_general_settings()
        with ui.CollapsableFrame("Branch Settings", name="group"):
            self.build_plant_branch_settings()
        with ui.CollapsableFrame("Leaf Settings", name="group"):
            self.build_plant_leaf_settings()
        with ui.HStack(spacing=HSPACING):
            ui.Label("Leaf Tex Name", name="attribute_name", width=LABEL_WIDTH)
            ui.StringField(self.plantSettings.leafTexture)
        with ui.HStack(spacing=HSPACING):
            ui.Label("Stem Tex Name", name="attribute_name", width=LABEL_WIDTH)
            ui.StringField(self.plantSettings.branchTexture)

    # ext_id is current extension id. It can be used with extension manager to query additional information, like where
    # this extension is located on filesystem.
    def on_startup(self, ext_id):
        print("[drillgon.scene.lsystemplants] MyExtension startup")

        self._window = ui.Window("Plant Generator", width=500, height=500)
        with self._window.frame:
            with ui.ScrollingFrame():
                with ui.VStack(height=0):
                    def on_line_draw():
                        self.isDrawingLine = True
                        self.linePointAnchor = None
                    
                    def on_plant_draw():
                        self.isDrawingPlant = True
                    
                    with ui.CollapsableFrame("Cables", name="group"):
                        with ui.VStack(height=0, spacing=VSPACING):
                            with ui.HStack(spacing=HSPACING):
                                myCurveTooltip = "My curves generate UVs differently, which may be desirable"
                                ui.Label("Use Built In Curves", name="attribute_name", width=LABEL_WIDTH, tooltip=myCurveTooltip)
                                ui.CheckBox(self.cableUseCurves, tooltip=myCurveTooltip)
                            with ui.HStack(spacing=HSPACING):
                                ui.Label("Curve Size", name="attribute_name", width=LABEL_WIDTH)
                                ui.FloatDrag(self.cableSize)
                            ui.Button("Draw Cable", clicked_fn=lambda: on_line_draw())
                    with ui.CollapsableFrame("Plants", name="group"):
                        with ui.VStack(height=0, spacing=VSPACING):

                            def reset_settings():
                                self.plantSettings.reset()

                            def reset_textures():
                                self.plantSettings.reset_textures()

                            self.build_plant_settings()
                            ui.Button("Reset Texture Names", clicked_fn=lambda: reset_textures())
                            ui.Button("Reset All Settings", clicked_fn=lambda: reset_settings())
                            with ui.HStack(spacing=HSPACING):
                                ui.Label("Plant Group", name="attribute_name", width=LABEL_WIDTH)
                                ui.StringField(self.plantSettings.group)

                            ui.Button("Draw Plant", clicked_fn=lambda: on_plant_draw())
                            


        appWindow = omni.appwindow.get_default_app_window()
        self._input = carb.input.acquire_input_interface()
        mouse = appWindow.get_mouse()
        keyboard = appWindow.get_keyboard()

        self._mouseEventHandler = self._input.subscribe_to_input_events(self.mouse_event_handler, device=mouse, order=0)
        self._keyEventHandler = self._input.subscribe_to_input_events(self.key_event_handler, device=keyboard, order=0)

    def on_shutdown(self):
        self._input.unsubscribe_to_input_events(self._mouseEventHandler)
        self._input.unsubscribe_to_input_events(self._keyEventHandler)
        print("[drillgon.scene.lsystemplants] MyExtension shutdown")
