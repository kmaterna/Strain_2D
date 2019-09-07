#!/bin/bash

range=-125/-120/38/42.1
projection=M5i
output="combo_vel_map.ps"

gmt pscoast -R$range -J$projection -B1.0 -Dh -Wthin,black -Sgray -Gwhite -N1 -N2 -K > $output


gmt psxy -R$range -J$projection -Wthick,black -K -O -P<<EOF>>$output
-123.260666 39.701538
#-126.167770 40.093243
#-125.803238 41.828274
-124.37 39.86
-123.98 41.6
-122.816429 41.437393
-123.260666 39.701538
EOF

gmt psxy ../../Vel_Data/Mapping_Resources/Quaternary.txt -R$range -J$projection -Wthinner,black -K -O >> $output


awk '{print $1, $2, $3, $4, $5, $6, $7}' CMM4_2006.gmtvec | gmt psvelo -R$range -J$projection -Se0.03/0.95/16 -Wred -Gred -K -O >> $output
awk '{print $1, $2, $3, $4, $5, $6, $7}' CEA_062806.gmtvec | gmt psvelo -R$range -J$projection -Se0.03/0.95/16 -Wblue -Gblue -K -O >> $output
awk '{print $2, $3, -10*$4, 10*$5, 10*$6, 10*$7, $8}' Williams_2006_dataset.txt | gmt psvelo -R$range -J$projection -Se0.03/0.95/16 -Worange -Gorange -K -O >> $output
awk '{print $1, $2, $3, $4, $5, $6, $7}' unr_velo.txt | gmt psvelo -R$range -J$projection -Se0.03/0.95/16 -Wpurple -Gpurple -K -O >> $output

gmt psvelo -R$range -J$projection -Se0.03/0.95/10 -Wred -Gred -K -O <<EOF >> $output
-124.4 39 20.0 0 0.5 0.5 0.0 CMM4
EOF
gmt psvelo -R$range -J$projection -Se0.03/0.95/10 -Wblue -Gblue -K -O <<EOF >> $output
-124.4 38.85 20.0 0 0.5 0.5 0.0 CEA
EOF
gmt psvelo -R$range -J$projection -Se0.03/0.95/10 -Worange -Gorange -K -O <<EOF >> $output
-124.4 38.7 20.0 0 0.5 0.5 0.0 Williams
EOF
gmt psvelo -R$range -J$projection -Se0.03/0.95/10 -Wpurple -Gpurple -K -O <<EOF >> $output
-124.4 38.55 20.0 0 0.5 0.5 0.0 UNR
EOF



rm gmt.history
open $output