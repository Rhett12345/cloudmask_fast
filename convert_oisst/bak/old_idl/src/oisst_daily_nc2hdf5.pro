pro oisst_daily_nc2hdf5

;---------------------------------------------------
; The propose of this program is to convert 
; oisst daily netcdf file to hdf5 format
;
; Author: Min Min
; Uint  : National Satellite Meteorological Center
;
;---------------------------------------------------

;--------------------------------------
;               start  
;--------------------------------------
print,'--------------------------------------'
print,'                start   '
print,'--------------------------------------'
print,' '
print,' Convert oisst daily netcdf file to hdf5 format !!! '

; Step [1] : read input filename using command line
print,' ' 
print,' Step [1] : read input information using command line '
args = command_line_args(count=count)
if count lt 3 then begin
   print, 'ERROR: args is wrong!'
   return
endif
nc_path   = args[0]
out_path  = args[1]
time1     = args[2]
year0  = strmid(time1,0,4)
month0 = strmid(time1,4,2)
day0   = strmid(time1,6,2)
print,'nc_path  = ',nc_path
print,'out_path = ',out_path
jday0 = julday(month0, day0, year0)
print,'satellite obs time = ', year0, month0, day0,' julian day =',jday0
jday = jday0 - 2
caldat, jday, month, day, year
print,'oisst obs time = ', year, month, day,' julian day =',jday
jday1 = julday(1, 1, year)
numday = jday - jday1 + 1
print,'number of day in one year =', numday

; Step [2] : find nc file name and define output hdf5 file name 
print,' ' 
print,' Step [2] : find nc file name and define output hdf5 file name '
nc_file = nc_path+'sst.day.mean.'+strtrim(string(year),2)+'.v2.nc'
print,'nc_file = ',nc_file
if (month lt 10) then begin
   mm = '0'+strtrim(string(month),2)
endif else begin
   mm = strtrim(string(month),2)
endelse
if (day lt 10) then begin
   dd = '0'+strtrim(string(day),2)
endif else begin
   dd = strtrim(string(day),2)
endelse
out_file = out_path+'sst.day.mean.'+strtrim(string(year),2)+mm+dd+'.hdf5'
print,'out_file = ',out_file


; Step [3] : read nc file 
print,' ' 
print,' Step [3] : read nc file  '
;sst0 = fltarr(1440,720,366)
nc_read, nc_file,'lon',sst0
n = size(sst0)
print,n

;--------------------------------------
;               end  
;--------------------------------------
print,' '
print,'--------------------------------------'
print,'                end   '
print,'--------------------------------------'

end


;--------------------------------------
;              Function s
;--------------------------------------

pro nc_read, filename, datasetName, data
    on_error, 2
    fid = ncdf_open(filename)
    varid = ncdf_varid(fid, datasetName)
    ncdf_varget, fid, varid, data
    ncdf_close, fid
end

pro hdf5_write, filename, datasetName, data, groupName=groupName, create=create
    if keyword_set(create) then begin
        if file_test(filename) then file_delete, filename
        fid = h5f_create(filename)
        if keyword_set(groupName) then groupID = h5g_create(fid, groupName)
    endif else begin
        if ~file_test(filename) then fid = h5f_create(filename) else $
            fid = h5f_open(filename, /write)
        if keyword_set(groupName) then begin
            nums = h5g_get_num_objs(fid)
            if nums eq 0 then groupID = h5g_create(fid, groupName) else begin
                for i=0, nums-1 do begin
                    name = h5g_get_obj_name_by_idx(fid, i)
                    if name eq groupName then break
                endfor
                if i eq nums then groupID = h5g_create(fid, groupName) else $
                    groupID = h5g_open(fid, groupName)
            endelse
        endif
    endelse
    
    datatypeID = h5t_idl_create(data)
    dataspaceID = h5s_create_simple(size(data, /dimensions))
    if keyword_set(groupName) then datasetID = h5d_create(groupID, datasetName, $
        datatypeId, dataspaceID) else $
        datasetID = h5d_create(fid, datasetName, datatypeId, dataspaceID)
    h5d_write, datasetID,data
    
    h5d_close, datasetID
    h5s_close, dataspaceID
    h5t_close, datatypeID
    if keyword_set(groupName) then h5g_close, groupID
    h5f_close, fid
end