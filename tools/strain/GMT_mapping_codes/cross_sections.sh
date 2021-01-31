#!/bin/bash
# Intended to be called from the Results/midas/Means folder

file1="../Results_hammond/I2nd.nc"
file2="../Results_spline/I2nd.nc"
file3="../Results_ND_interp/I2nd.nc"
file4="../Results_visr/I2nd.nc"
file5="../Results_tape/I2nd.nc"
file6="../Results_gpsgridder/I2nd.nc"
file7="means_I2nd.nc"
tempfile="tempfile.xyz"
output1="projection1.txt"
output2="projection2.txt"
output3="projection3.txt"
output4="projection4.txt"
output5="projection5.txt"
output6="projection6.txt"
output7="projection7.txt"
range="-125/-121.2/39.0/42.0"
cx=-122.55
cy=40
azimuth=75
l_min=-110
l_max=110
l_span=220
w_min=-3
w_max=3
w_span=6
projection="M3.8i"
plot_range="-110/110/-1/5"
plot_size="4.5i/2.8i"
out_ps="I2nd_coplot.ps"
out_pdf="I2nd_coplot.pdf"


rm $out_pdf

# projecting data to line
gmt grd2xyz $file1 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output1
rm $tempfile

gmt grd2xyz $file2 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output2
rm $tempfile

gmt grd2xyz $file3 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output3
rm $tempfile

gmt grd2xyz $file4 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output4
rm $tempfile

gmt grd2xyz $file5 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output5
rm $tempfile

gmt grd2xyz $file6 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output6

gmt grd2xyz $file7 > $tempfile
gmt project $tempfile -C$cx/$cy -L$l_min/$l_max -A$azimuth -W$w_min/$w_max -Q -Fpz > $output7

# plotting california
gmt makecpt -T-1/5/0.1 -D -Crainbow.cpt > cpt1.cpt
gmt psbasemap -R$range -J$projection -Y6 -BWESN+t"Cross-Section Location" -Bp1.0 -K > $out_ps
gmt grdimage means_I2nd.nc -R$range -J$projection -BWeSN -Bp1.0 -Ccpt1.cpt -K -O >> $out_ps
gmt pscoast -R$range -J$projection -Wthick,black -Df -Slightblue -K -O >> $out_ps

# adding location of xsection
gmt psxy -R$range -J$projection -SJ -K -O -Wthickest <<EOF>> $out_ps
$cx $cy $azimuth $l_span $w_span
EOF

# # plotting profiles of strain
gmt psbasemap -R$plot_range -JX$plot_size -Y3 -X13 -Bxf10a50 -Byf1a2 -Bx+l"Km" -By+l"nanostrain(log)" -BwESn+t"Second Invariant" -K -O >> $out_ps
gmt psxy $output1 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,red3 -Gred3 -K -O >> $out_ps
gmt psxy $output2 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,indianred1 -Gindianred1 -K -O >> $out_ps
gmt psxy $output3 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,goldenrod1 -Ggoldenrod1 -K -O >> $out_ps
gmt psxy $output4 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,limegreen -Glimegreen -K -O >> $out_ps
gmt psxy $output5 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,royalblue -Groyalblue -K -O >> $out_ps
gmt psxy $output6 -R$plot_range -JX$plot_size -Sc0.04i -W0.1p,purple3 -Gpurple3 -K -O >> $out_ps
# adding lines to graph
# sort -r -k1,1 -n $output1 | gmt psxy -R$plot_range -JX$plot_size -W0.1p,red3 -K -O >> $out_ps
# sort -r -k1,1 -n $output2 | gmt psxy -R$plot_range -JX$plot_size -W0.1p,indianred1 -K -O >> $out_ps
# sort -r -k1,1 -n $output3 | gmt psxy -R$plot_range -JX$plot_size -W0.1p,goldenrod1 -K -O >> $out_ps
# sort -r -k1,1 -n $output4 | gmt psxy -R$plot_range -JX$plot_size -W0.1p,limegreen -K -O >> $out_ps
# sort -r -k1,1 -n $output5 | gmt psxy -R$plot_range -JX$plot_size -W0.1p,royalblue -K -O >> $out_ps
# # plotting means
sort -r -k1,1 -n $output7 | gmt psxy -R$plot_range -JX$plot_size -Wthick,black -K -O >> $out_ps

gmt pslegend -F+gazure1+pblack -Dx0.i/-2.2i+w2.1i/1.6i+jBL+l1.2 -O -K <<EOF>> $out_ps
G 0.05i
N 1
S 0.1i c 0.1i red3 0.25p 0.3i Delaunay
S 0.1i c 0.1i indianred1 0.25p 0.3i Numpy Spline
S 0.1i c 0.1i goldenrod1 0.25p 0.3i Numpy ND Interp
S 0.1i c 0.1i limegreen 0.25p 0.3i VISR
S 0.1i c 0.1i royalblue 0.25p 0.3i Tape/Wavelets
S 0.1i c 0.1i purple3 0.25p 0.3i GPSgridder
S 0.1i - 0.1i 0.25i black 0.3i Mean
EOF
rm $output1
rm $output2
rm $output3
rm $output4
rm $output5
rm $output6
rm $output7
rm gmt.history
gmt psconvert $out_ps -Tf