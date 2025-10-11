load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;;;************************************************
begin

data=asciiread("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/mg_1950-2022.txt",-1,"integer")
;;;;各类路径气旋的位置
npts=dimsizes(data)
print("Number of values: "+npts)
x=data
dim = dimsizes(x)

f  = addfile("/cygdrive/e/data/track_data/TRACK_ff_trs_pos_ERA5_GRID1_1950101_20230720_NH_CYC.nc","r")
;;;;;气旋数据
parcel=f->parcel
time=f->traj_time 
vor=f->traj_rel_vor_pos
time_FillValue=-2147483647
time1=time*6
time1@units="hours since 1950-01-01 00:00:0.0"
time1@long_name="Time of the trajectory"
Time=cd_calendar(time1(:,0),-5)

vor@_FillValue=-999.0 

x1=ind(Time(x,0).ge.1950.and.Time(x,0).le.2022)
x1=x(x1)
print(dimsizes(x1))
;;;MAM
x2=ind(Time(x,1).ge.3.and.Time(x,1).le.5)
x2=x(x2)
nn2=dimsizes(x2)
print(nn2)

;;;JJA
x3=ind(Time(x,1).ge.6.and.Time(x,1).le.8)
x3=x(x3)
nn3=dimsizes(x3)
print(nn3)

;;;SON
x4=ind(Time(x,1).ge.9.and.Time(x,1).le.11)
x4=x(x4)
nn4=dimsizes(x4)
print(nn4)

;;;DJF
x5=ind(Time(x,1).ge.12.or.Time(x,1).le.2)
x5=x(x5)
nn5=dimsizes(x5)
print(nn5)

;;;;;
XM=x5
a=dimsizes(XM)

;;;;;;;统计每条气旋路径轨迹个数
numb=new(a,float)
numb@_FillValue=-9.999
pvor=new(a,float)
pvor@_FillValue=-9.999
;;;;;vor
do i = 0,a-1
numb(i)=num(.not.ismissing(vor(XM(i),:)))
 pvor(i)=max(vor(XM(i),:));;;;;;max
 ; pvor(i)=vor(x(i),0);;;;;;start
 ; pvor(i)=avg(vor(x(i),:));;;;;;avg
end do
midvor=dim_median(pvor);;;;中位数
; ;;;;;;;;;;maxvor>3
 ; p=ind(pvor.ge.midvor)
 ; p=x(p)

; ;;;;;按照强度排序
 ; nn=floattoint(a/3.0);;;;;;;;;;;;;
; xppvor=new((/a,2/),float)
; ipp=dim_pqsort(pvor,1);;;;;;;;;;;1递增、-1递减
; xppvor(:,0)=x(ipp)
; xppvor(:,1)=pvor(ipp)
; xpp_90=xppvor(0:nn,0)

; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/mg_maxvor_-33.txt",xpp_90)
; exit
delete(vor)
;;;;;每年生成的气旋各自对应的长度与速度
lon=f->traj_lon
lat=f->traj_lat
lon@_FillValue=-999.0
lat@_FillValue=-999.0 

lat_sta=lat(XM,:)
lon_sta=lon(XM,:)
len=new((/a,150/),float)
len@_FillValue=-9.96921e+36
; len(:,0)=0.0
length=new(a,float)
length@_FillValue=-9.96921e+36
sp=new((/a,150/),float)
sp@_FillValue=-9.96921e+36
; sp(:,0)=0.0
speed=new(a,float)
speed@_FillValue=-9.96921e+36
do i = 0,a-1
do j = 1,149
if(.not.ismissing(lat_sta(i,j)))
len(i,j)=gc_latlon(lat_sta(i,j-1),lon_sta(i,j-1), lat_sta(i,j),lon_sta(i,j), 2,4)
length(i)=sum(len(i,:));;;km
sp(i,j) = len(i,j)*1000/(6*3600);;;m/s
end if
end do
end do
  do i=0,a-1
  speed(i)=max(sp(i,:))
 end do

;;;;;;;;;;;size;;
delete(lat)
delete(lon)
size=f->traj_size
msize=new(a,float) 
msize@_FillValue=-2147483647

dd=new((/a,150/),"float")
   do i=0,a-1
     do j=0,149
dd(i,j)=sqrt(size(XM(i),j)*6378*6378/3.14)
	 end do
   end do  
  do i=0,a-1
   msize(i)=max(dd(i,:))
    end do

; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJF_mg_maxsize.txt",q)
; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJF_mg_life.txt",numb)
; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJF_mg_length.txt",length)
; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJF_mg_maxspeed.txt",speed)
; asciiwrite("/cygdrive/E/data/ERA5_TRACK/grid1/1950-2022/DJF_mg_maxvor.txt",pvor)

end