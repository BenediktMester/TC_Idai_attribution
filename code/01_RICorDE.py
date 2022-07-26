'''
Created on Mar. 26, 2022

@author: cefect


ricorde calcs for TC Idai data from AER

'''

import os, datetime, copy
start =  datetime.datetime.now()
print('start at %s'%start)
 


from ricorde.scripts import Session, QgsCoordinateReferenceSystem, force_open_dir

def runr(
        tag = 'r1',
        name = 'idai',
        crsid = 'EPSG:32737',
        aoi_fp= r'C:\LS\02_WORK\NRC\2202_TC\04_CALC\aoi\aoi04_0326.gpkg',
        dem_fp = r'C:\LS\10_OUT\2202_TC\ins\dem\merit_0304\MERIT_merge_0304_90x90_aoi04.tif',            
        inun_fp = r'C:\LS\02_WORK\NRC\2202_TC\06_DATA\aer\220307\aer_afed_hilo_acc_3s_20190301-20190331_v05r01_0326_xfed.tif',
        pwb_fp = r'C:\LS\02_WORK\NRC\2202_TC\06_DATA\JRC\JRC_extent_merge_0326_aoi05_clean.tif', #native resolution
 
        compress='med',
        
        #run_dataPrep
        pwb_resampling='Maximum',
        
        #build_b1Bounds: hand value stats for bouding beach1 samples
       qhigh=0.8, cap=6.0,  #uppers               
       qlow=0.2, floor=1.0, #lowers
        
        #build_inun1
        buff_dist=0, #pwb has lots of noise
        
       
        
        #build_beach2
        b2_method='polygons', b2_spacing=90*4, b2_write=True,
        
        #build_hgInterp
        hgi_minPoints=3, searchRad=90*12, hgi_resolution=90*6,
        
        #build_hgSmooth
        hval_precision=0.5,   max_iter=5,
        
        #build_depths
        d_compress='med', 
        
        **kwargs):
    
 
    
    with Session(name=name, tag=tag,
                 root_dir=r'C:\LS\10_OUT\2202_TC',
                 compress=compress,  
                 crs=QgsCoordinateReferenceSystem(crsid),
                   overwrite=True,
                   bk_lib = {
                       'pwb_rlay':dict(resampling=pwb_resampling),
                       'b1Bounds':dict(qhigh=qhigh, cap=cap, qlow=qlow, floor=floor),
                       'inun1':dict(buff_dist=buff_dist),
                       'beach2':dict(method=b2_method, spacing=b2_spacing, write_plotData=b2_write),
                       'hgInterp':dict(pts_cnt=hgi_minPoints, radius=searchRad, resolution=hgi_resolution),
                       'hgSmooth':dict(max_iter=max_iter, precision=hval_precision),
                       'depths':dict(compress=d_compress),

                       },
                   aoi_fp=aoi_fp, dem_fp=dem_fp, inun_fp=inun_fp, pwb_fp=pwb_fp, #filepaths for this project
                   **kwargs) as wrkr:
        
 

        #=======================================================================
        # wrkr.run_dataPrep()
        # wrkr.run_HAND()
        # wrkr.run_imax()
        # wrkr.run_HANDgrid()
        # wrkr.run_wslRoll()
        #=======================================================================
        wrkr.run_depths()
        
        out_dir = wrkr.out_dir
        
    return out_dir


def dev_run(
        tag = 'dev',
        floor=1.5,cap=7.0,
        
        **kwargs):
 
    
    return runr(
        tag=tag,floor=floor, cap=cap,
        
        aoi_fp= r'C:\LS\02_WORK\NRC\2202_TC\04_CALC\aoi\t\aoiT02_0327.gpkg', 
        dem_fp = r'C:\LS\10_OUT\2202_TC\ins\dem\merit_0304\MERIT_merge_0304_90x90_aoiT01.tif', 
        inun_fp = r'C:\LS\02_WORK\NRC\2202_TC\06_DATA\aer\220307\aer_afed_hilo_acc_3s_20190301-20190331_v05r01_0326_xfed_aoiT01.tif',
        pwb_fp = r'C:\LS\02_WORK\NRC\2202_TC\06_DATA\JRC\JRC_extent_merge_0326_aoiT02_clean.tif',

  
        b2_write=True,
        compress='none', 
        
        pwb_resampling='Maximum',
        
        qhigh=0.8, qlow=0.2, 
        max_iter=3, hval_precision=0.2,
        
        compiled_fp_d = {
        #=======================================================================
        # 'dem':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_dem.tif',
        # 'pwb_rlay':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_pwb_rlay.tif',
        # 'inun_rlay':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_inun_rlay.tif',
        # 'dem_hyd':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_dem_hyd.tif',
        # 'HAND':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\idai_dev_0331_HAND.tif',
        # 'HAND_mask':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_HAND_mask.tif',
        # 'inun1':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_inun1.tif',
        # 'beach1':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_beach1.tif',
        # 'b1Bounds':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_b1Bounds.pickle',
        # 'inunHmax':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_inunHmax.tif',
        # 'inun2':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\idai_dev_0331_inun2.tif',
        # 'beach2':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_beach2.gpkg',
        # 'hgInterp':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_hgInterp.tif',
        # 'hgRaw':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_hgRaw.tif',
        # 'hgSmooth':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\idai_dev_0331_hgSmooth.tif',
        # 'hInunSet':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_hInunSet.pickle',
        # 'hWslSet':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\working\idai_dev_0331_hWslSet.pickle',
        # 'wslMosaic':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\idai_dev_0331_wslMosaic.tif',
        #=======================================================================
        #'depths':r'C:\LS\10_OUT\2202_TC\outs\dev\20220331\idai_dev_0331_depths.tif',

            },
        
        **kwargs)


def r1(**kwargs):
    return runr(
        tag='r1',
#===============================================================================
#             compiled_fp_d={
#                 'pwb_rlay':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_pwb_rlay.tif',
#                 'inun_rlay':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_inun_rlay.tif',
#                 'dem':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_dem.tif',
#                 
#                 'dem_hyd':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_dem_hyd.tif',
#                 'HAND':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\idai_r1_0329_HAND.tif',
#                 'HAND_mask':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_HAND_mask.tif',
#                 
#                 'inun1':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_inun1.tif',
#                 'beach1':r'C:\LS\10_OUT\2202_TC\outs\r1\20220329\working\idai_r1_0329_beach1.tif',
#                 
#                 'b1Bounds':r'C:\LS\10_OUT\2202_TC\outs\r1\20220330\working\idai_r1_0330_b1Bounds.pickle',
#                 'inunHmax':r'C:\LS\10_OUT\2202_TC\outs\r1\20220330\working\idai_r1_0330_inunHmax.tif',
#                 'inun2':r'C:\LS\10_OUT\2202_TC\outs\r1\20220330\idai_r1_0330_inun2.tif',
#                 'beach2':r'C:\LS\10_OUT\2202_TC\outs\r1\20220330\working\idai_r1_0330_beach2.gpkg',
#                 
#                 'wslMosaic':r'C:\LS\10_OUT\2202_TC\outs\r1\20220331\idai_r1_0331_wslMosaic.tif',
# 
#                 },
#===============================================================================
        **kwargs)

if __name__ =="__main__": 
    
 
    output = dev_run()
    #output= r1()
    
 
    
    
    #===========================================================================
    # wrap
    #===========================================================================
    #force_open_dir(od)
    tdelta = datetime.datetime.now() - start
    print('finished in %s\n    %s'%(tdelta, output))
