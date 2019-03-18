try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

RecAA_4x7261ERHIO200NMSWAMP_rec_vtk = LegacyVTKReader( FileNames=['/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/261_263_273/Rec-AA_4x7-261-ERHIO200-NM-SW/Rec-AA_4x7-261-ERHIO200-NM-SW-AMP.rec.vtk'] )

RenderView1 = GetRenderView()
DataRepresentation1 = Show()
DataRepresentation1.BackfaceDiffuseColor = [0.40000000000000002, 0.40000000000000002, 1.0]
DataRepresentation1.Representation = 'Outline'
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation1.DiffuseColor = [0.40000000000000002, 0.40000000000000002, 1.0]

RenderView1.CameraPosition = [-7.9857177734375, 1.4235916137695312, 2126.2456493404034]
RenderView1.CameraClippingRange = [1511.9956171634055, 2912.8715795639077]
RenderView1.CameraFocalPoint = [-7.9857177734375, 1.4235916137695312, -3.663116455078125]
RenderView1.CameraParallelScale = 551.260952918675
RenderView1.CenterOfRotation = [-7.9857177734375, 1.4235916137695312, -3.663116455078125]

Contour1 = Contour( PointMergeMethod="Uniform Binning" )

Contour1.PointMergeMethod = "Uniform Binning"
Contour1.ContourBy = ['POINTS', 'scalars']
Contour1.Isosurfaces = [0.50001390183974037]

DataRepresentation2 = Show()
DataRepresentation2.BackfaceDiffuseColor = [0.40000000000000002, 0.40000000000000002, 1.0]
DataRepresentation2.Texture = []
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation2.DiffuseColor = [0.40000000000000002, 0.40000000000000002, 1.0]

RenderView1.CameraViewUp = [0.0, 0.0, 1.0]
RenderView1.CameraPosition = [669.36931048501049, 0.018096923828125, 7.7463760375976562]
RenderView1.CameraClippingRange = [461.16451134683518, 961.06114977149707]
RenderView1.CameraFocalPoint = [-12.378498077392578, 0.018096923828125, 7.7463760375976562]
RenderView1.CameraParallelScale = 176.44931681285726
RenderView1.CenterOfRotation = [-12.378498077392578, 0.018096923828125, 7.7463760375976562]

DataRepresentation1.Visibility = 0

DataRepresentation2.DiffuseColor = [0.0, 1.0, 0.50196078431372548]

WriteImage('/Users/jesseclark/Documents/MATLAB/data_analysis/Au0710/261_263_273/Rec-AA_4x7-261-ERHIO200-NM-SW/test.jpg')


Render()
