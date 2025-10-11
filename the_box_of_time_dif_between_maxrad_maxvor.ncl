;
;  These files still have to be loaded manually
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/box_percentile_plot.ncl"
;********************************************
begin 
data1=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/MAMmg_1950-2022.txt",-1,"integer")
 ; data1=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/MAMmg_1950-2022.txt",-1,"float")
 x1=data1

 data2=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/JJAmg_1950-2022.txt",-1,"integer")
 x2=data2

 data3=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/SONmg_1950-2022.txt",-1,"integer")
 x3=data3

 data4=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJFmg_1950-2022.txt",-1,"integer")
 x4=data4

 f  = addfile("/cygdrive/e/data/track_data/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")
;气旋数据
time=f->traj_time 
vor=f->traj_rel_vor_pos
vor@_FillValue=-999.0 
size=f->traj_size
time_FillValue=-2147483647

y=x1
nnumb1=dimsizes(y)
;;;;;;判断涡度最大的时间点
p=new(nnumb1,float) 
p@_FillValue=-999
do i=0,nnumb1-1
   p(i)=maxind(vor(y(i),:))
 end do
;;;;;;;;;判断半径最大的时间点
q=new(nnumb1,float) 
q@_FillValue=-2147483647

dd=new((/nnumb1,150/),"float")
   do i=0,nnumb1-1
     do j=0,149
dd(i,j)=sqrt(size(y(i),j)*6378*6378/3.14)
	 end do
   end do
do i=0,nnumb1-1

   q(i)=maxind(dd(i,:))
 end do

 ;;;;;;;;;
 m=new(nnumb1,float) 
m@_FillValue=-999.0
m=p-q

 
dimm1=dimsizes(data1)
dimm2=dimsizes(data2)
dimm3=dimsizes(data3)
dimm4=dimsizes(data4)

dat                           = new((/4,dimm1/),float)
   dat@_FillValue                =   32767
   
   dat(0,:)                      =  m ;	dif_lon_avg;   
  delete(p)  
    delete(q)  
	 delete(dd) 
	  delete(m) 
    delete(y)  
  delete(nnumb1) 
  
y=x2
nnumb1=dimsizes(y)
;;;;;;判断涡度最大的时间点
p=new(nnumb1,float) 
p@_FillValue=-999
do i=0,nnumb1-1
   p(i)=maxind(vor(y(i),:))
 end do
;;;;;;;;;判断半径最大的时间点
q=new(nnumb1,float) 
q@_FillValue=-2147483647

dd=new((/nnumb1,150/),"float")
   do i=0,nnumb1-1
     do j=0,149
dd(i,j)=sqrt(size(y(i),j)*6378*6378/3.14)
	 end do
   end do
do i=0,nnumb1-1

   q(i)=maxind(dd(i,:))
 end do

 ;;;;;;;;;
 m=new(nnumb1,float) 
m@_FillValue=-999.0
m=p-q
   dat(1,0:dimm2-1)                      =  m ;	dif_lon_avg; 
  delete(p)  
    delete(q)  
	 delete(dd) 
	  delete(m) 
    delete(y)  
  delete(nnumb1) 
   
  y=x3
nnumb1=dimsizes(y)
;;;;;;判断涡度最大的时间点
p=new(nnumb1,float) 
p@_FillValue=-999
do i=0,nnumb1-1
   p(i)=maxind(vor(y(i),:))
 end do
;;;;;;;;;判断半径最大的时间点
q=new(nnumb1,float) 
q@_FillValue=-2147483647

dd=new((/nnumb1,150/),"float")
   do i=0,nnumb1-1
     do j=0,149
dd(i,j)=sqrt(size(y(i),j)*6378*6378/3.14)
	 end do
   end do
do i=0,nnumb1-1

   q(i)=maxind(dd(i,:))
 end do

 ;;;;;;;;;
 m=new(nnumb1,float) 
m@_FillValue=-999.0
m=p-q
   dat(2,0:dimm3-1)                      =  m ;	dif_lon_avg; 
  delete(p)  
    delete(q)  
	 delete(dd) 
	  delete(m) 
    delete(y)  
  delete(nnumb1) 
  
  y=x4
nnumb1=dimsizes(y)
;;;;;;判断涡度最大的时间点
p=new(nnumb1,float) 
p@_FillValue=-999
do i=0,nnumb1-1
   p(i)=maxind(vor(y(i),:))
 end do
;;;;;;;;;判断半径最大的时间点
q=new(nnumb1,float) 
q@_FillValue=-2147483647

dd=new((/nnumb1,150/),"float")
   do i=0,nnumb1-1
     do j=0,149
dd(i,j)=sqrt(size(y(i),j)*6378*6378/3.14)
	 end do
   end do
do i=0,nnumb1-1

   q(i)=maxind(dd(i,:))
 end do

 ;;;;;;;;;
 m=new(nnumb1,float) 
m@_FillValue=-999.0
m=p-q
   dat(3,0:dimm4-1)                      =  m ;	dif_lon_avg; 
  delete(p)  
    delete(q)  
	 delete(dd) 
	  delete(m) 
    delete(y)  
  delete(nnumb1) 
     
	
  mmin=min(dat)
  mmax=max(dat)
	
  wks = gsn_open_wks("eps","/cygdrive/E/program4/1the_box_of_time_dif_maxrad_maxvor")

  res                          = True

  res@tmXBLabels               = (/"MAM","JJA","SON","DJF"/);"JRA55";;"clus1","clus2","clus3","clus4";;"point0","p_maxvor"
  res@tmXBLabelFontHeightF     = 0.02 
  res@tmXBLabelAngleF          = 0
  res@tmXBLabelJust            = "CenterCenter"
  res@tmYMajorGrid                = True;   True       ; implement y grid
  res@tmYMajorGridLineDashPattern = True             ; select short dash lines
  res@trYMinF                  =    -20;mmin-10; 
  res@trYMaxF                  =   20;mmax+10;2000;;1500;;50; ; 
  res@trXMinF                  = -1 ;-90
  res@trXMaxF                  = 4 ; 90
  res@tiMainFontHeightF        = 0.023
  ; res@tiMainString             = "Example Box-Percentile-Plot"
  res@tmXTBorderOn             = True
  res@tmYRBorderOn             = True 
  res@tmYUseLeft              = True
  res@vpXF                     = 0.15
  res@vpYF                     = 0.9
  res@vpWidthF                 = 0.85
  res@vpHeightF                = 0.65 
  
  boxOpts                      = True
  boxOpts@boxWidth             = 0.3 ; with DNA == True 0.3 / without 0.4
  boxOpts@boxColors          = "red"
  lineRes                      = True
  lineRes@gsLineThicknessF     = 0.5
  lineRes@gsLineColor          ="blue"
  lineRes@txFontHeightF            = 5 
  lineRes@tmXBLabelFontHeightF     = 9
  markerRes                    = True 
  markerRes@gsMarkerIndex      = 1
  markerRes@gsMarkerSizeF      = 0.05
  markerRes@gsMarkerThicknessF = 1
  markerRes@gsMarkerColor      = "black"; (/"","blue"/) ;"red"
  PerRes                       = True
  PerRes@gsLineColor           = "black"; () ;
  PerRes@gsLineThicknessF      = 2  
  PerRes@gsFillColor           = (/"red3","seagreen3","royalblue3","orange"/) ; "grey";rot
  PerRes@FillBox               =  True;(/True,True,/);
  PerRes@MarkPer               = True ;(/True,True/)
  PerRes@Avg                   = True ;(/True,False/)
  DNA                          = False; True;;;;;例图中的右侧多个小横线
  
  NoMember                     = True ; False ;True;;;;;;;;;;;每列个数
  NoMember@txFontHeightF       = 0.015 
  NoMember@txJust              = "BottomLeft"
  box   = box_percentile_plot( wks,dat, res,boxOpts,lineRes,PerRes,markerRes,DNA,NoMember)
  txres = True
  gsres = True
  txres@txFontHeightF = 0.015 
  txres@txJust = "CenterCenter"
  
  txres@txJust = "CenterCenter"
  txres@txFontHeightF = 0.02 
  gsn_text_ndc(wks,"point",0.1,0.200,txres) 

  draw(box)
   frame(wks)
 ;************************************************
 end ; program
;************************************************   
	