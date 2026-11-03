load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function attach_labelbar(wks,map,arr[*]:numeric,colors[*])
local lbres, vph, vpw, nboxes
begin
  getvalues map
    "vpHeightF" : vph
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)
   
  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbFillColors       = colors
  lbres@vpWidthF           = 0.7 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.2 * vph     ; labelbar height
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.01          ; font height. default is small
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False
  lbres@lbLabelAlignment   = "InteriorEdges"
  lbid = gsn_create_labelbar(wks,nboxes,""+(arr),lbres)

  amres                  = True
  amres@amJust           = "TopCenter"
  amres@amParallelPosF   =  0.0    ; Center
  amres@amOrthogonalPosF =  0.6    ; Bottom
  annoid = gsn_add_annotation(map,lbid,amres)
  return(annoid)
end

begin

cy=addfile("/cygdrive/E/data/weather/wp.1954-2015.nc","r")

data=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/MAMmg_1950-2022.txt",-1,"integer")
data1=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/JJAmg_1950-2022.txt",-1,"integer")
data2=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/SONmg_1950-2022.txt",-1,"integer")
data3=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJFmg_1950-2022.txt",-1,"integer")

f=addfile("/cygdrive/E/data/ERA5_TRACK/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")

;;;;各类路径气旋的位置
npts=dimsizes(data)
print("Number of values: "+npts)
x=data
dim1= dimsizes(x)
   delete(data)
path=ispan(0,dim1-1,1)
point=ispan(0,149,1)

;;;;;;;;;;气旋数据
time=f->traj_time 
tlon=f->traj_lon
tlat=f->traj_lat
se=f->traj_size
time_FillValue=-2147483647
t1=time*6
t1@units="hours since 1950-01-01 00:00:0.0"
t1@long_name="Time of the trajectory"
Time1=cd_calendar(t1(x,0:149),-2)
tlon_p=tlon(x,:)
tlat_p=tlat(x,:)
size=se(x,:)
anumb=new(dim1,integer)
anumb@_FillValue=-2147483647 
do i = 0,dim1-1
anumb(i)=num(.not.ismissing(tlon_p(i,:)))
end do
dd=new((/dim1,150/),"float")
dd@_FillValue=-999.0 
   do i=0,dim1-1
     do j=0,149
dd(i,j)=sqrt(size(i,j)*6378*6378/3.14)
	 end do
   end do

; print(Time1)
   delete(tlat)
   delete(tlon)
   delete(t1)
   delete(f)
   delete(se)
   delete(time)
   delete(size)
   
   
;;;;;;;;;;;;;;;;;;降水数据
sn=cy->station
sd=cy->dss
sr=cy->rd
sf=cy->fd
lon=cy->lon
lat=cy->lat
t2=cy->time
sd@_FillValue=8
sr@_FillValue=8
sf@_FillValue=8
ssm=sd+sr+sf
ssm@_FillValue=8
print(max(ssm))
ssm=where(ssm.ge.1,1,ssm);;;;;;;三种天气现象相加，出现一种以上计为1
ssm@_FillValue=8
copy_VarMeta(sd,ssm)


time_units = "hours since 1900-1-1 00:00:0.0"
t2@units = time_units
Tim2=cd_calendar(t2,-5)
yt=ind(Tim2(:,0).ge.1955.and.Tim2(:,0).le.2014)
ss=ssm(:,yt)

; print(Tim2)
delete(ssm)
delete(Tim2)
delete(cy)
printVarSummary(ss)
 
tt=ss&time
tt@units = time_units
Time2=cd_calendar(tt,-2)
    delete(tt)

dim2=dimsizes(sn)
dim3=dimsizes(Time2)
day=ispan(0,34,1)
printVarSummary(t2)
; exit
ssd=new((/dim3,dim2/),float)
ssd=int2flt(ss(time|:,station|:))
              printVarSummary(ssd)
delete(ss)


hs=new(dim2,float)			
do i=0,dim2-1
hs(i)=num(ssd(:,i).eq.1)
end do  


 
;;;;;;;;;;;;;;;;气旋大小范围内站点的累积降水量

aveds=new((/dim1,150,dim2/),integer)
aveds@_FillValue=8
ddsn=new(dim2,float)
 do i=0,dim1-1
   do j=0,149;;;;;气旋的时间
    ddsn   =  gc_latlon(tlat_p(i,j),tlon_p(i,j), lat(:),lon(:),2,4)
    num2   =   num(dd(i,j).ge.ddsn)
   if(num2.gt.0)
    id2  =  ind(dd(i,j).ge.ddsn)
   
    aveds(i,j,id2)  = (/sn(id2) /)
    delete(id2)

  end if
    delete(num2)
    delete(ddsn)
 end do
end do	
   delete(tlat_p)
   delete(tlon_p)
   
   
dsssn=new((/dim3,dim2/),float)
dsssn@_FillValue=8
  do i=0,dim1-1
 ntime	=	get_unique_values(Time1(i,:))
 nnt=dimsizes(ntime)
 
  do k=0,nnt-1

	idx1	=	ind(Time1(i,:).eq.ntime(k));;;;将同一日的轨迹点挑出来	
	nnst    =   get_unique_values(aveds(i,idx1,:));;;;统计某一天气旋轨迹内的总站点
	num1    =   num(nnst)
	num2	=	num(ntime(k).eq.Time2);;;;对应站点数据的时间
	
   if(num1.gt.0.and.(num2.gt.0))
	id_dif  =   venn2_difference(floattoint(nnst),sn);;;;挑出轨迹内站点对应到所有站点的下标
    idx2    =   get1Dindex_Exclude(sn,id_dif)
	idx3	=	ind(ntime(k).eq.Time2);;;;对应站点数据的时间
	; avedss(i,k,idx2)   =   ss(idx2,idx3)
    dsssn(idx3,idx2)   =   ssd(idx3,idx2);+dsssn(idx3,idx2)
; print(dsssn(idx3,idx2))
; exit
   delete(id_dif)
   delete(idx2)
   delete(idx3)
   end if
   delete(idx1)
   delete(num1)
   delete(num2)
   delete(nnst)
  end do  
   delete(ntime)
   delete(nnt)
end do  
copy_VarMeta(ssd(:,:),dsssn(:,:))

 
printVarSummary(dsssn)
outfil1=addfile("/cygdrive/E/data/wp_MAMmg_1955-2014.nc","c") 
outfil1->sss=dsssn
 outfil1@history="sd+sr+sf DJF mg 1955-2014"
 exit
;;;;;;;CDO:monthly sum
;;;;;;;;;;season sum
 dssMonM=new((/62,dim2/),float)
 dssMonM@_FillValue=8
do i=0,61
dssMonM(i,:)=dim_sum_n_Wrap(dssMon(12*i+2:12*i+4,:),0)
end do
dssMonMAM=dim_avg_n_Wrap(dssMonM,0)
dssMonMAM@_FillValue=8
 outf1=addfile("/cygdrive/E/data/wp_clus1_sumMAM.nc","c") 
 outf1->wp=dssMonMAM
 outf1@history="wp sumMAM clus1 1954-2015 avg"

; exit

;;;;;;;;;;;;挑出没有产生沙尘暴的气旋
; indss=new(dim1,float)
; indss@_FillValue=32766
; indss=dim_sum_n_Wrap(dim_sum_n_Wrap(avedss,2),1)
; x1=ind(.not.ismissing(indss))
; x1=x(x1)
; printVarSummary(x1)
 ; asciiwrite("/public/home/wangji19/data/allmeg/spr/have_wp_es.txt",x1)
  
; ssdss=dim_avg_n_Wrap(dsssn,0)
 ; ssdss=month_to_season(dsssn,"MAM")  
; ssdss=dim_sum_n_Wrap(dim_sum_n_Wrap(avedss(:,:,:),1),0)
 
print(max(dssMonMAM))
print(min(dssMonMAM))
 
; ; ssdss=where(ssdss.eq.0,32766,ssdss)
; dssMAM@_FillValue=32766

; exit
; print(ssdss)
;;;;;;;;;;;;;;;;;;;;气旋生命期内所有站点的降水
; preci=new((/dim1,dim2/),float)
; preci@_FillValue=-999.0
; do i = 0,dim1-1;;;;气旋个数
; id2=ind(Time1(i).eq.Time2)
; if(.not.ismissing(id2)) 
; preci(i,:)=dim_sum_n_Wrap(pre(:,id2:anumb(i)+id2-1),1)
; end if
; end do 

; preci=where(ismissing(preci),0,preci)
; spre=new(dim2,float)
; spre@_FillValue=-999.0
; spre=dim_avg_n_Wrap(preci(:,:),0)


;;;;;;;;;;;;;;;;;;;70年所有气旋累积
arr=(/1,2.,3.,4.,5.,6.,13./)
; colors = (/4,5,6,7,8,9,10,11/)
colors = (/168,180,193,206,220,230,242,255/)

;;;;;;;;;;;;;;;;;;;;70年所有气旋平均
; arr=(/0.5,1.,1.5,2.,2.5,3./)
; colors = (/49,81,113,145,177,209,255/)

num_markers=dimsizes(arr)+1;
; colors = (/2,6,10,14,18,22,26,30,34,38,42,46,50,54,60/)
;colors  = span_color_rgba("MPL_PuBu",num_markers)
;colors  = span_color_rgba("MPL_RdBu",num_markers)
lat_new = new((/num_markers,dimsizes(dssMonMAM)/),float);
lon_new = new((/num_markers,dimsizes(dssMonMAM)/),float);

;;;;;---------------------------
do i =0,num_markers-1
        if(i.eq.0);第一个水平
            indexes=ind(dssMonMAM.lt.arr(0))
        end if
        if(i.eq.(num_markers-1))then;最大的一个水平
            indexes=ind(dssMonMAM.ge.max(arr))
        end if
        if(i.gt.0.and.i.lt.(num_markers-1))then;中间的水平
            indexes=ind(dssMonMAM.ge.arr(i-1).and.dssMonMAM.lt.arr(i))
        end if
        ;把lat、lon_new的前N个（在这一水平里有N个点）换成这N个点的经纬度。
    if (.not.any(ismissing(indexes))) then
            npts_range=dimsizes(indexes)
            lat_new(i,0:(npts_range-1))=lat(indexes)
            lon_new(i,0:(npts_range-1))=lon(indexes)
    end if

        delete(indexes)
end do

; print(lat_new)
; exit
wks=gsn_open_wks("eps","/cygdrive/E/data/station_sumMAM_avg_wp_of_the_clus1(1955-2014)")
gsn_define_colormap(wks, "BlueWhiteOrangeRed");GMT_drywetsunshine_9lev

  res                            = True            
  res@gsnMaximize                = False
  res@gsnDraw                    = False
  res@gsnFrame                   = False
  res@pmTickMarkDisplayMode      = "Always"
  ; res@gsnLeftString="wp during the life of es path(1954-2014)"
  res@gsnLeftStringOrthogonalPosF=0.02
  ; res@tiXAxisString="mm"
  res@tiXAxisFontHeightF=0.017
  res@tiXAxisOffsetYF=-0.07

;---------------------------------------------------------------
; set for the map 
  res@mpMinLatF                  = 17.
  res@mpMaxLatF                  = 55.
  res@mpMinLonF                  = 72.
  res@mpMaxLonF                  = 136.
   
  res@mpFillOn                   = True
  res@mpDataSetName              = "/cygdrive/E/data/MAP/NCL-Chinamap-master/database/Earth..4"
  res@mpDataBaseVersion          = "MediumRes" ; or "Ncarg4_1"
  res@mpAreaMaskingOn            = True
  res@mpMaskAreaSpecifiers       = (/"China"/)
  res@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
    
  res@mpLandFillColor            = "white"
  res@mpInlandWaterFillColor     = "white"
  res@mpOceanFillColor           = "white"
  ;res@mpOutlineBoundarySets      = "NoBoundaries"
  res@mpNationalLineColor        = "black"
  res@mpProvincialLineColor      = "black"
  res@mpGeophysicalLineColor     = "black"
  res@mpNationalLineThicknessF   = 2
  res@mpProvincialLineThicknessF = 1

  map = gsn_csm_map(wks,res)  
  
  res1=True
  res1@gsLineColor="black"
  plot1=gsn_add_shapefile_polylines(wks,map,filename,res1)
    
 
  gsres               = True
  gsres@gsMarkerIndex = 1         ; Use filled dots for markers.

  ;base_size = 0.01
  pmid = new(num_markers,graphic)
  do i = 0, num_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i)
	  ;gsres@gsMarkerColor      = colors(i,:)
      gsres@gsMarkerSizeF      = 0.015;base_size * (i+1)/3;
      ;gsres@gsMarkerThicknessF = 1.5*(i+1)
      dumstr = unique_string("marker")
      map@$dumstr$ = gsn_add_polymarker(wks,map,lon_new(i,:),lat_new(i,:),gsres)
	  
	  end if

  end do
  lbid =   attach_labelbar(wks,map,arr,colors)
  ; draw(map)
  ; frame(wks)
  
  
;--- add South China Sea --- 
 nhres                          = res
  nhres@gsnMaximize              = False

  nhres@vpHeightF                = 0.13    
  nhres@vpWidthF                 = 0.13
    
  nhres@mpMinLatF                =   2.0    
  nhres@mpMaxLatF                =  20.0
  nhres@mpMinLonF                = 107.0
  nhres@mpMaxLonF                = 123.0

  ;nhres@lbLabelBarOn             = False
  nhres@tmXBOn                   = False 
  nhres@tmXTOn                   = False
  nhres@tmYLOn                   = False
  nhres@tmYROn                   = False
  nhres@gsnLeftString            = ""
  nhres@gsnRightString           = ""
  nhres@tiXAxisString=""

  map_nanhai = gsn_csm_map(wks,nhres)
  plot2=gsn_add_shapefile_polylines(wks,map_nanhai,filename,res1)
  pmid1 = new(num_markers,graphic)
  do i = 0, num_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i)
   ; gsres@gsMarkerColor      = colors(i,:)
      gsres@gsMarkerSizeF      = 0.015;base_size * (i+1)/3;
      ; gsres@gsMarkerThicknessF = 1.5*(i+1)
      pmid1(i) = gsn_add_polymarker(wks,map_nanhai,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

  adres                          = True
  adres@amParallelPosF           = 0.47 ; -0.5 is the left edge of the plot.
  adres@amOrthogonalPosF         = 0.47  ; -0.5 is the top edge of the plot.
  adres@amJust                   = "BottomRight"

  plotnh = gsn_add_annotation(map,map_nanhai,adres)
  draw(map)
  frame(wks)
  
  
end

