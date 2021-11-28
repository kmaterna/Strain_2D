#!/bin/bash
# GMT plotting of several strain results

range=$1
projection="M2.0i"
out_strain="output_rows.ps"
gmt set MAP_FRAME_TYPE plain
gmt set FORMAT_GEO_MAP D

output_dir="../Output/"
infile1=$output_dir"delaunay/delaunay_strain.nc"
infile2=$output_dir"gpsgridder/gpsgridder_strain.nc"
infile3=$output_dir"visr/visr_strain.nc"
infile4=$output_dir"loc_avg_grad/loc_avg_grad_strain.nc"

# Delaunay Strain
gmt makecpt -T-1/5/0.5 -Cbatlow.cpt > mycpt.cpt
gmt grdedit $infile1=gd?HDF5:"$infile1"://I2 -R$range -T -G$"I2_delaunay.nc" 
gmt grdedit I2_delaunay.nc -Ev
gmt grdimage I2_delaunay.nc -R$range -J$projection -BWeSn -Bp1.0 -Cmycpt.cpt -X2 -Y10 -K > $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
echo "-122.8 42.3 Delaunay" | gmt pstext -R -J -F+f18p,Helvetica-Bold -N -K -O >> $out_strain
gmt psvelo $output_dir"delaunay/positive_eigs.txt" -Se0.003/0.68/0 -A+e+n10+pthick,blue -Gblue -R -J -K -O >> $out_strain
gmt psvelo $output_dir"delaunay/negative_eigs.txt" -Se0.003/0.68/0 -A+b+n10+pthick,black -Gred -R -J -K -O >> $out_strain

# Loc_Avg_Grad Strain
gmt grdedit $infile4=gd?HDF5:"$infile4"://I2 -R$range -T -G$"I2_loc_avg_grad.nc" 
gmt grdedit I2_loc_avg_grad.nc -Ev
gmt grdimage I2_loc_avg_grad.nc -R$range -J$projection -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
echo "-122.8 42.3 Loc avg grad" | gmt pstext -R -J -F+f18p,Helvetica-Bold -N -K -O >> $out_strain
gmt psvelo $output_dir"loc_avg_grad/positive_eigs.txt" -Se0.003/0.68/0 -A+e+n10+pthick,blue -Gblue -R -J -K -O >> $out_strain
gmt psvelo $output_dir"loc_avg_grad/negative_eigs.txt" -Se0.003/0.68/0 -A+b+n10+pthick,black -Gred -R -J -K -O >> $out_strain

# Visr Strain
gmt grdedit $infile3=gd?HDF5:"$infile3"://I2 -R$range -T -G$"I2_visr.nc" 
gmt grdedit I2_visr.nc -Ev
gmt grdimage I2_visr.nc -R -J -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
echo "-122.8 42.3 VISR" | gmt pstext -R -J -F+f18p,Helvetica-Bold -N -K -O >> $out_strain
gmt psvelo $output_dir"visr/positive_eigs.txt" -Se0.003/0.68/0 -A+e+n10+pthick,blue -Gblue -R -J -K -O >> $out_strain
gmt psvelo $output_dir"visr/negative_eigs.txt" -Se0.003/0.68/0 -A+b+n10+pthick,black -Gred -R -J -K -O >> $out_strain

# GPSgridder Strain
gmt grdedit $infile2=gd?HDF5:"$infile2"://I2 -R$range -T -G$"I2_gpsgridder.nc" 
gmt grdedit I2_gpsgridder.nc -Ev
gmt grdimage I2_gpsgridder.nc -R -J -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
echo "-122.8 42.3 gpsgridder" | gmt pstext -R -J -F+f18p,Helvetica-Bold -N -K -O >> $out_strain
gmt psscale -DjTR+w4.5/0.5+o-1.1/1.5 -R -J -B1.0:"log(I2)":/:: -Cmycpt.cpt -K -O >> $out_strain
gmt psvelo $output_dir"gpsgridder/positive_eigs.txt" -Se0.003/0.68/0 -A+e+n10+pthick,blue -Gblue -R -J -K -O >> $out_strain
gmt psvelo $output_dir"gpsgridder/negative_eigs.txt" -Se0.003/0.68/0 -A+b+n10+pthick,black -Gred -R -J -K -O >> $out_strain


# Rotation
gmt makecpt -T0/300/10 -Cmagma.cpt -Do -G0.30/1.0 > mycpt.cpt
gmt grdedit $infile1=gd?HDF5:"$infile1"://rotation -R$range -T -G$"rot_delaunay.nc" 
gmt grdedit rot_delaunay.nc -Ev
gmt grdimage rot_delaunay.nc -R -J -BWeSn -Bp1.0 -Cmycpt.cpt -X-18 -Y-8 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
awk '{print $1, $2, $3, $4, $6, $7, 0}' $output_dir"delaunay/tempgps.txt" | gmt psvelo -Se0.03/0.68/0 -A+e+pthick,black -Gwhite -R -J -K -O >> $out_strain

gmt grdedit $infile4=gd?HDF5:"$infile4"://rotation -R$range -T -G$"rot_loc_avg_grad.nc" 
gmt grdedit rot_loc_avg_grad.nc -Ev
gmt grdimage rot_loc_avg_grad.nc -R -J -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
awk '{print $1, $2, $3, $4, $6, $7, 0}' $output_dir"loc_avg_grad/tempgps.txt" | gmt psvelo -Se0.03/0.68/0 -A+e+pthick,black -Gwhite -R -J -K -O >> $out_strain

gmt grdedit $infile3=gd?HDF5:"$infile3"://rotation -R$range -T -G$"rot_visr.nc" 
gmt grdedit rot_visr.nc -Ev
gmt grdimage rot_visr.nc -R -J -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
awk '{print $1, $2, $3, $4, $6, $7, 0}' $output_dir"visr/tempgps.txt" | gmt psvelo -Se0.03/0.68/0 -A+e+pthick,black -Gwhite -R -J -K -O >> $out_strain

gmt grdedit $infile2=gd?HDF5:"$infile2"://rotation -R$range -T -G$"rot_gpsgridder.nc" 
gmt grdedit rot_gpsgridder.nc -Ev
gmt grdimage rot_gpsgridder.nc -R -J -BweSn -Bp1.0 -Cmycpt.cpt -X6 -Y0 -K -O >> $out_strain
gmt pscoast -R -J -Wthick,black -Df -Sgray -K -O >> $out_strain
awk '{print $1, $2, $3, $4, $6, $7, 0}' $output_dir"gpsgridder/tempgps.txt" | gmt psvelo -Se0.03/0.68/0 -A+e+pthick,black -Gwhite -R -J -K -O >> $out_strain
gmt psscale -DjTR+w4.5/0.5+o-1.1/1.5 -R$range -J$projection -B50:"Rotation":/:: -Cmycpt.cpt -O >> $out_strain

rm gmt.history
rm gmt.conf
rm mycpt.cpt
gmt psconvert -Tg $out_strain

open $out_strain
