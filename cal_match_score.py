# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:39:22 2016

@author: kezhili
"""
import math
import numpy as np
import scipy.io as sio
import sys

def cal_match_score(cvs_ind, mov_fra_ind, mov_fra_ind_beg, para_extend, csv_mag, frame_mag):
    # initialize
    leng_cvs_ind = cvs_ind.shape[0]
    leng_csv = csv_mag.shape[0]
    weights = np.power(math.e,np.arange(0,-2.0001,-(2/(leng_cvs_ind-1))))
    
    # match_mtx is used to record all possible penalties for different mm1 and
    # mm2
    match_mtx = np.zeros((leng_cvs_ind,4),dtype=float);
    match_mtx[:,0] = cvs_ind
    
    mov_fra_ind_extend = np.round_((mov_fra_ind[(mov_fra_ind_beg-1):]-mov_fra_ind[mov_fra_ind_beg-1])*para_extend + cvs_ind[0])
    
    # parameters to calculate the 'shift_to_left'
    shift_weights = np.power(math.e,np.arange(-2.4,0.01,0.6))
    shift_weights = shift_weights.conj().transpose()
    normal_shift_weights = shift_weights/shift_weights.sum(axis=0)
    # only conder nearest 10 peaks to calculate the 'shift_to_left'
    shift_consider = np.zeros((1,shift_weights.shape[0]),dtype=int)
    
    
    # compensate the difference: shift impulse index to left 
    shift_to_left = np.zeros((leng_cvs_ind+1,1),dtype=int)
    
    
    
    mov_fra_ind_extend_copy = mov_fra_ind_extend.copy()
    
    for ii in range(1,leng_cvs_ind+1):
        dis_seq1 =  abs(mov_fra_ind_extend_copy-cvs_ind[ii-1]-shift_to_left[ii-1])
        min_ind = np.argmin(dis_seq1) +1       
        
        match_ind = mov_fra_ind_extend_copy[min_ind-1];
        match_mtx[ii-1,1] = match_ind;
        
        # update 'shift_to_left'. it is caluclated based on last shift values
        shift_to_left[ii-1] = match_mtx[ii-1,1] -match_mtx[ii-1,0];
        shift_consider[0:-1] = shift_consider[1:];
        #shift_consider(end) = min_fra_ind_match(ii,5)-mask_ind_x(ii,1);
        shift_consider[-1] =  shift_to_left[ii-1];
        shift_to_left_val =np.dot(shift_consider,normal_shift_weights)
        shift_to_left[ii] = round(shift_to_left_val[0])
        
        mov_fra_ind_extend_copy[min_ind-1] = 0    
    
#    for ii = 1:leng_cvs_ind;
#        match_mtx(ii,3) = csv_mag(cvs_ind(ii));
#        match_mtx(ii,4) = sum(frame_mag(max(1,match_mtx(ii,2)-3):min(match_mtx(ii,2)+3,leng_csv)));
#    end
    for ii in range(1,leng_cvs_ind+1):
        match_mtx[ii-1,2] = csv_mag[cvs_ind[ii-1]-1]
        surr_ind = np.arange( max(1,match_mtx[ii-1,1]-3), min(match_mtx[ii-1,1]+3,leng_csv)+1,dtype = int)
        surr_val = frame_mag[surr_ind-1];        
        match_mtx[ii-1,3] = surr_val.sum(axis=0)
    
    zero_ind = np.where(match_mtx[:,1]==0)[0]+1
    positive_ind = np.where(match_mtx[:,1]>0)[0]+1
    
    #if sum(abs(sort(match_mtx(:,2))-match_mtx(:,2)))>1
    if min(zero_ind, default= 1e6) < max(positive_ind, default=1):
        match_score = 1e10
    else:
        # absulute difference/ and/ second derivitive 
        match_diff = match_mtx[:,1] - match_mtx[:,0]
        
         # distance from mm1=1, mm2=20
        mm1_add = mov_fra_ind_beg-1
        mm2_add = abs(round((para_extend-0.94)/0.003)-20)
        
        # calculate the match score, given weights, and distance mm1=1, mm2=20
        weights1 = abs(match_diff*weights)   
        weights1 = weights1.sum(axis=0)
        weights2 = mm1_add/2+ mm2_add*2
        match_score = weights1 * (1 + weights2*0.02)
        
    return match_score, match_mtx

if __name__ == "__main__":
#    cvs_ind = np.array([  274, 1749,3356,6185,6274,6362,6488,6872,6999,7513,8193,8867,9898,9994, 10309,10522,10702,11945,12264,13628,
#                        14524,14764,14901,15096,15188,15313, 15468,15550,15622,16045,16127,16219,16696,16735,16770,16828,16927,17007,
#                        17068,17109, 17154,17199,17241, 17280,17321,17366,17412,17456,17498,17594,17898, 18030,18229,19370,19591,19796,
#                        19921,19993,20045,20094,20157,20222,20410,20521, 20586, 20797,20834,20874,20923,20973,21018, 21061,21107,21147,
#                        21213,21302, 21408,21606,21712,21809,22071,22405, 22453,22511,22622,22684,23038,23101,23199,23246,23283,
#                        23388,23416,23444,23484,23526, 23558,23598,23651,23699,23748, 23792, 23837,23882,23959,24098,24173,24235,
#                        24295,24585,24687,24763,24937,25058,25317,25414,25950,26073,26443, 26747])
#                        
#    mov_fra_ind = np.array([ 273,1746,3350,6178,6267,6355,6480, 6864,6991, 7505,8183, 8857,9886, 9982,10298,10509, 10689,11932,
#                             12250,13613, 14507, 14747,14884, 15079, 15170,15296,15450,15533,15604,16027, 16108, 16201, 16677,
#                             16716,16751,16808,16907, 16988, 17049,17089,17134,17179, 17222,17261,17301,17346,17392, 17436,17479,
#                             17574,17878, 18010,18207,19347, 19569,19774, 19898,19972,20022, 20071, 20135, 20199,20386,20498,
#                             20563, 20773, 20811, 20850,20899,20948,20994, 21037,21083,21123,21188,21277,21384,21582,21688, 21784,
#                             22046,22380, 22427,22486,22596, 22658,23012, 23074,23172, 23219,23256,23360,23389, 23418, 23457,
#                             23499, 23532, 23571, 23624, 23672,23721,23765,23810, 23855, 23932,24072,24146,24207,24267,24557,
#                             24659,24735,24908,25030,25289,25385,25920,26044, 26413,26717])
#    mov_fra_ind_beg = 1
#    para_extend =1
#    diff_mask_central=sio.loadmat('diff_mask_central_temp1.mat')
#    diff_mask_central = diff_mask_central['diff_mask_central']
    
    cal_match_read = sio.loadmat('cal_match_read.mat')
    cvs_ind = cal_match_read['cvs_ind']
    mov_fra_ind = cal_match_read['mov_fra_ind']
    mov_fra_ind_beg = cal_match_read['mm1']
    para_extend = cal_match_read['mm2']*0.003+0.94
    diff_mask_central = cal_match_read['diff_mask_central']
    
    csv_mag=  diff_mask_central[:,6]                  
    frame_mag = diff_mask_central[:,2]                    
    
    match_score, match_mtx  = cal_match_score(cvs_ind.flatten(), mov_fra_ind.flatten(), mov_fra_ind_beg[0][0], para_extend[0][0], csv_mag, frame_mag)
#    match_score = cal_match_score(cvs_ind, mov_fra_ind, mov_fra_ind_beg, para_extend, csv_mag, frame_mag)
    
    sys.stdout.write(str(match_score)) 
    sys.stdout.write(str(match_mtx))
    
    sio.savemat('match_mtx_res.mat',mdict = {'match_score':match_score, 'match_mtx':match_mtx})
    
#    sio.savemat('cal_match_read.mat',mdict = {'cvs_ind':cvs_ind, 'mov_fra_ind':mov_fra_ind, 'mm1': 1, 'mm2': 20, 'diff_mask_central': diff_mask_central})