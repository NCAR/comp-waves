
;PURPOSE: File to store manual date information
;         Also contains default options
;
;
;          Modularised by R J Morton 2015
;          Created temp_index to keep together useful quantities from wave tracking RJM 04/2016
;
;
PRO find_man_date,date ;date,x1,x2,y1,y2,start_file,$
                   ;num_files_proc,lower_r,upper_r,ccthresh,maxiter

COMMON wt_var,inpath,nwlst,rsun_mm, temp_index,ccthresh,maxiter,num_files_proc

; *** Default values for a date that has not ***
; *** been checked before and is not defined ***
; *** below in the case statement            ***

; coordinates for cross-correlation box
temp_index.cc_coord=[[45,230],[75,390]] ;[x1,y1],[x2,y2]


; Choose start file - used to skip bad files
temp_index.start_file = 0

; Number of files to process. 0 means: let the 
; code choose the highest number possible. Useful 
; for manually cutting out bad data at the end.
num_files_proc = 0

; Upper and lower limit for the radius
; !be aware that the maximum offset of the cross-correlation
; (if performed) is added to this value!
temp_index.lower_r = 228.
temp_index.upper_r = 288.

; threshold for the cross-correlation precision
ccthresh = 0.04

; maximum iteration steps for the cross-correlation
maxiter = 200. 


CASE date OF
 ;example for a manual setting for a specific date

 '20111230': begin
                temp_index.cc_coord=[[540,319],[575,371]]
                temp_index.start_file = 0
                num_files_proc = 179
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
               end

 '20120327': begin
               temp_index.cc_coord=[[55,220],[85,260]]
      	       temp_index.start_file = 0
      	       num_files_proc = 164
      	       temp_index.lower_r = 228.
      	       temp_index.upper_r = 280.
              end

 '20120411': begin
              temp_index.cc_coord=[[540,350],[580,390]]
      	      temp_index.start_file = 0        ;361
      	      num_files_proc = 360  ;532
      	      temp_index.lower_r = 228.
      	      temp_index.upper_r = 290.
             end

 '20120410': begin
               temp_index.cc_coord=[[535,390],[565,430]]
      	       temp_index.start_file = 0
      	       num_files_proc = 300
      	       temp_index.lower_r = 228.
      	       temp_index.upper_r = 290.
      	      end

 '20120706': begin
               temp_index.cc_coord=[[539,282],[563,310]]
               temp_index.start_file = 42
               num_files_proc = 132
               temp_index.lower_r = 228.
               temp_index.upper_r = 290.
             end

 '20130213': BEGIN
              temp_index.cc_coord=[[545,375],[580,410]]
              temp_index.start_file = 0
              num_files_proc = 178
              temp_index.lower_r = 228.
              temp_index.upper_r = 288.
             END

 '20130502': begin
                temp_index.cc_coord=[[540,220],[570,250]]
      	        temp_index.start_file = 0
      	        num_files_proc = 176
      	        temp_index.lower_r = 228.
      	        temp_index.upper_r = 290.
            end

 '20130515': begin
                temp_index.cc_coord=[[545,331],[573,372]]
                temp_index.start_file = 0
                num_files_proc = 173
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end

 '20130708': begin
                temp_index.cc_coord=[[57,241],[80,275]]
                temp_index.start_file = 0
                num_files_proc = 179
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end                      

 '20130914': begin
                temp_index.cc_coord=[[545,310],[580,345]]
                temp_index.start_file = 0
      	        num_files_proc = 152
      	        temp_index.lower_r = 228.
      	        temp_index.upper_r = 290.
             end

 '20131223': begin
                temp_index.cc_coord=[[61,178],[92,223]]
                temp_index.start_file = 0
                num_files_proc = 48
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end

 '20140511': begin
                temp_index.cc_coord=[[61,199],[89,231]]
                temp_index.start_file = 0
                num_files_proc = 178
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
              end

 '20140510': begin
                temp_index.cc_coord=[[538,345],[570,377]]
                temp_index.start_file = 0
                num_files_proc = 180
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end

  '20150121': begin
                temp_index.cc_coord=[[46,220],[77,259]]
                temp_index.start_file = 0
                num_files_proc = 169
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end
  '20150122': begin
                temp_index.cc_coord=[[46,220],[77,259]]
                temp_index.start_file = 0
                num_files_proc = 169
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end                       

  '20150202': begin
                temp_index.cc_coord=[[77,140],[112,176]]
                temp_index.start_file = 0
                num_files_proc = 169
                temp_index.lower_r = 228.
                temp_index.upper_r = 290.
             end

  else:      begin
               message, /cont, 'Manual settings for this date were not found, using default values.'
             end
endcase



END